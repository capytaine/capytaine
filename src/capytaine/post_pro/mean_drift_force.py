# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Computation of the mean drift force."""

import logging

import numpy as np
import xarray as xr

from capytaine.bem.airy_waves import (
    airy_waves_pressure,
    airy_waves_velocity,
)
from capytaine.bodies.bodies import FloatingBody
from capytaine.bodies.multibodies import Multibody
from capytaine.bodies.dofs import AbstractDof
from capytaine.bem.problems_and_results import DiffractionProblem
from capytaine.tools.perturbation_expansion import DifferenceFrequencyExpansion as Expansion


LOG = logging.getLogger(__name__)


####################################################################################################
# Far field
####################################################################################################

def compute_far_field_mean_drift_force(dataset, *, rao):
    """Compute the horizontal mean drift force using the far field formulation.

    The force is computed from the Kochin functions of the diffracted and radiated waves.
    Only the horizontal components (surge, sway and yaw) are available.

    Parameters
    ----------
    dataset : xarray Dataset
        The dataset of the diffraction and radiation results, including the Kochin functions
        ``kochin_diffraction`` and ``kochin_radiation`` (computed when ``theta`` is in the test matrix).
        The range of ``theta`` should cover [0, 2pi], and strictly include the wave directions.
    rao : xarray DataArray
        The complex motion amplitude of the body per unit wave amplitude,
        as returned by :func:`capytaine.post_pro.rao.rao`.
        The frequency should be its first dimension.

    Returns
    -------
    xarray Dataset
        A dataset with a ``mean_drift_force_far_field`` variable of dimensions
        ``(freq, wave_direction, other_wave_direction, influenced_dof)``,
        where ``freq`` is the frequency dimension of ``rao`` and ``influenced_dof`` is
        ``["Surge", "Sway", "Yaw"]``.
        It follows the same convention as :func:`compute_near_field_mean_drift_force`:
        it is Hermitian in the two wave directions, and the term of the pair of waves
        ``(wave_direction=k, other_wave_direction=l)`` is multiplied by ``a_k conj(a_l)``,
        where ``a_k`` is the complex amplitude of the wave ``k``.
    """
    omega = dataset['omega']
    m = dataset['wavenumber']
    h = dataset['water_depth']
    beta = dataset['wave_direction']
    rho = dataset['rho']
    g = dataset['g']
    H_diff = dataset['kochin_diffraction']
    H_rad = dataset['kochin_radiation']
    theta_range = dataset['theta']

    if (theta_range.min() > 0) or (theta_range.max() < 2*np.pi):
        raise ValueError("Theta should takes values between 0 and 2pi")
    if np.any(beta < theta_range.min()) or np.any(beta > theta_range.max()):
        raise ValueError("The wave direction should be in the theta interval")
    if np.any(beta == theta_range.min()) or np.any(beta == theta_range.max()):
        raise ValueError("The wave direction should not be at the border of the theta interval, it is recommended to extend the theta interval")

    H_rad_tot = sum(H_rad.sel(radiating_dof=d)*rao.sel(radiating_dof=d) for d in rao.radiating_dof)
    H_tot = np.exp(1j*np.pi/2)*(H_diff + H_rad_tot)

    H_beta = H_tot.interp(theta=beta.values)
    H_derivative = H_tot.differentiate("theta")
    H_derivative_beta = H_derivative.interp(theta=beta.values)

    coef1 = 2*np.pi*rho*omega
    if h == np.inf:
        coef2 = 2*np.pi*rho*m**2
    else:
        k0 = omega**2/g
        coef2 = 2*np.pi*rho*m*(k0*h)**2 / (h*((m*h)**2 - (k0*h)**2 + k0*h))

    dims = [rao.dims[0], "wave_direction", "other_wave_direction"]
    coords = {
        rao.dims[0]: rao.coords[rao.dims[0]].values,
        "wave_direction": rao.coords["wave_direction"].values,
        "other_wave_direction": rao.coords["wave_direction"].values,
    }

    freq = rao.sizes[rao.dims[0]]
    nb_dir = rao.sizes["wave_direction"]
    data_x = np.full((freq, nb_dir, nb_dir), np.nan + 1j*np.nan, dtype=complex)
    data_y = np.full((freq, nb_dir, nb_dir), np.nan + 1j*np.nan, dtype=complex)
    data_z = np.full((freq, nb_dir, nb_dir), np.nan + 1j*np.nan, dtype=complex)
    H = H_tot.sel(theta=slice(0, 2*np.pi))
    H_derivative = H_derivative.sel(theta=slice(0, 2*np.pi))

    for k in range(nb_dir):
        for l in range(k, nb_dir):
            hk = H.sel(wave_direction=beta[k])
            hl = H.sel(wave_direction=beta[l])
            hk_derivative = H_derivative.sel(wave_direction=beta[k])
            hl_derivative = H_derivative.sel(wave_direction=beta[l])

            data_x[:, k, l] = (-coef1 * (np.cos(beta[l]) * H_beta.sel(theta=beta[l], wave_direction=beta[k]) - np.cos(beta[k]) * np.conjugate(H_beta.sel(theta=beta[k], wave_direction=beta[l]))) / (2*1j)
                               - coef2 * (hk* np.conjugate(hl) * np.cos(H.theta)).integrate("theta"))
            data_x[:, l, k] = np.conjugate(data_x[:, k, l])
            data_y[:, k, l] = (-coef1 * (np.sin(beta[l]) * H_beta.sel(theta=beta[l], wave_direction=beta[k]) - np.sin(beta[k]) * np.conjugate(H_beta.sel(theta=beta[k], wave_direction=beta[l]))) / (2*1j)
                               - coef2 * (hk* np.conjugate(hl) * np.sin(H.theta)).integrate("theta"))
            data_y[:, l, k] = np.conjugate(data_y[:, k, l])
            data_z[:, k, l] = (coef1/m * (H_derivative_beta.sel(theta=beta[l], wave_direction=beta[k]) + np.conjugate(H_derivative_beta.sel(theta=beta[k], wave_direction=beta[l]))) / 2
                               - coef2/m * ((np.conjugate(hk) * hl_derivative - hl * np.conjugate(hk_derivative)) / (2*1j)).integrate("theta"))
            data_z[:, l, k] = np.conjugate(data_z[:, k, l])

    force = xr.DataArray(
        data=np.stack([data_x, data_y, data_z], axis=-1),
        dims=[*dims, "influenced_dof"],
        coords={**coords, "influenced_dof": ["Surge", "Sway", "Yaw"]},
        )
    return xr.Dataset({"mean_drift_force_far_field": force})


_LEGACY_FAR_FIELD_VARIABLES = {"Surge": "drift_force_surge", "Sway": "drift_force_sway", "Yaw": "drift_force_yaw"}
_LEGACY_FAR_FIELD_DIMS = {"wave_direction": "wave_direction_k", "other_wave_direction": "wave_direction_l"}


def far_field_mean_drift_force(X, dataset):
    """Deprecated alias of :func:`compute_far_field_mean_drift_force`, with the output format of Capytaine 3.0.

    Parameters
    ----------
    X : xarray DataArray
        The motion RAO.
    dataset : xarray Dataset
        This function supposes that variables named 'kochin_diffraction' and 'kochin_radiation' are in the dataset.

    Returns
    -------
    xarray Dataset
        A dataset with the variables ``drift_force_surge``, ``drift_force_sway`` and ``drift_force_yaw``,
        of dimensions ``(freq, wave_direction_k, wave_direction_l)``.
    """
    LOG.warning("Deprecation warning: `far_field_mean_drift_force` is deprecated, "
                "please use `compute_far_field_mean_drift_force`, which returns a single variable "
                "`mean_drift_force_far_field` with the dimensions `wave_direction`, `other_wave_direction` and `influenced_dof`.")
    force = compute_far_field_mean_drift_force(dataset, rao=X)["mean_drift_force_far_field"]
    force = force.rename(_LEGACY_FAR_FIELD_DIMS)
    return xr.Dataset({name: force.sel(influenced_dof=dof, drop=True).rename(name)
                       for dof, name in _LEGACY_FAR_FIELD_VARIABLES.items()})


def _far_field_mean_drift_force_from_legacy_variables(dataset):
    """Inverse of the formatting of `far_field_mean_drift_force`: the variables `drift_force_surge`,
    `drift_force_sway` and `drift_force_yaw` merged into a DataArray in the format of
    `compute_far_field_mean_drift_force`, for datasets computed with Capytaine 3.0."""
    force = xr.concat(
        [dataset[name] for name in _LEGACY_FAR_FIELD_VARIABLES.values()],
        dim=xr.DataArray(list(_LEGACY_FAR_FIELD_VARIABLES), dims=["influenced_dof"], name="influenced_dof"),
    )
    force = force.transpose(..., "influenced_dof")
    return force.rename({v: k for k, v in _LEGACY_FAR_FIELD_DIMS.items()}).rename("mean_drift_force_far_field")


####################################################################################################
# Near field
####################################################################################################

def compute_near_field_mean_drift_force(
        dataset, body,
        *, rao, solver=None, output_pressure=False
):
    """Compute the mean drift force using the near field formulation.

    The force is the second-order part of the integral of the pressure on the moving body and along
    the waterline, for any dofs (rigid body motions or generalized modes with a gradient of motion).
    Only the interaction between waves at the same frequency (with the same or with different wave directions)
    are considered and all the terms of the result are at difference-frequency zero,
    where the second-order potential does not contribute to the pressure.

    The first-order fields are evaluated from a dataset computed with the indirect method and
    ``keep_details=True`` (e.g. ``solver.fill_dataset(test_matrix, body, keep_details=True)``),
    which stores the pressures and the source distributions on the body.
    Several frequencies and dofs are taken into account at once.
    Only datasets with a single value of rho, g and water depth, and without forward speed, are supported.

    Parameters
    ----------
    dataset : xarray Dataset
        The dataset of the diffraction and radiation results, for all the wave frequencies
        and wave directions of interest, with the source distributions on the body
        (``sources_diffraction``, ``sources_radiation``, and on the lid if any).
    body : FloatingBody or Multibody
        The body of the dataset.
    rao : xarray DataArray
        The complex motion amplitude of the body per unit wave amplitude,
        as returned by :func:`capytaine.post_pro.rao.rao`.
    solver : BEMSolver, optional
        The solver evaluating the potential and its gradient from the source distributions.
        It is recommended to use the same one as the one used to compute the dataset.
        Default: initialize default ``BEMSolver()``.
    output_pressure : bool, optional
        If True, also return the second-order pressure field on the hull and the first-order
        fields on the waterline from which the waterline term of the force is computed.
        Default: False.

    Returns
    -------
    xarray Dataset
        A dataset with a ``mean_drift_force_near_field`` variable of dimensions
        ``(freq, wave_direction, other_wave_direction, influenced_dof)``,
        where ``freq`` is the frequency dimension of the inputs (e.g. ``omega`` or ``wavelength``).
        It is Hermitian in the two wave directions: the term of the pair of waves
        ``(wave_direction=k, other_wave_direction=l)`` is multiplied by ``a_k conj(a_l)``,
        where ``a_k`` is the complex amplitude of the wave ``k``.
        If the inputs have different frequencies, wave direction or dofs, the output is
        computed on their intersection.
        If ``output_pressure`` is True, the dataset also contains the following variables:

        - ``pressure_hull_mean_drift``, of dimensions
          ``(freq, wave_direction, other_wave_direction, hull_face)``: the mean drift
          pressure field on the hull. It has the same normalization as
          ``mean_drift_force_near_field``: its diagonal is the mean pressure in a single wave of
          unit amplitude. It does not include the waterline term.
        - ``pressure_waterline_total``, of dimensions
          ``(freq, wave_direction, waterline_edge)``: the complex first-order
          pressure ``p`` in the strip between the mean waterline and the free surface, at the
          middle of each edge of the waterline, for a wave of unit amplitude.
          It is the total pressure on the freely floating body: incident, diffracted and
          radiated (with the motion given by ``rao``) waves,
          as well as the hydrostatic pressure averaged over the height of the strip.
        - ``relative_elevation_waterline``, of the same shape: the complex first-order relative
          elevation ``h`` of the free surface with respect to the body, at the same points, for
          a wave of unit amplitude.
          These two first-order fields depend on a single wave direction. Their mean product
          for a pair of waves can be computed by renaming ``wave_direction`` to
          ``other_wave_direction`` in a copy of them (see cookbook example B10).
        - ``mesh_waterline_vertices``, of shape ``(nb_waterline_edges, 2, 3)``: the coordinates of
          the two ends of each edge of the waterline.
    """
    if solver is None:
        from capytaine.bem.solver import BEMSolver  # Not at the top of the module, to avoid a circular import.
        solver = BEMSolver()
    body = body.as_FloatingBody if isinstance(body, Multibody) else body
    mesh = body.mesh.merged()  # Symmetries are ignored for now
    _check_inputs_for_near_field(dataset, body, solver)
    rao, dataset, body, freq_type = _align_coordinates(
            rao, dataset, body)
    motion_amplitude = rao.values  # (nb_freq, nb_wave_direction, nb_dofs)
    rho = float(dataset.coords["rho"])
    g = float(dataset.coords["g"])

    # Still water (order 0): hydrostatic pressure on the hull and geometry of the waterline
    p0 = -rho * g * mesh.faces_centers[:, 2]  # (nb_faces,)
    grad_p0 = np.array([0.0, 0.0, -rho * g])
    vertical_dofs_motion_waterline = body._dofs_on_waterline()[:, :, 2]  # (nb_dofs, nb_edges_waterline)
    inclination_waterline = np.sqrt(1 - mesh.waterline_edges_normals[:, 2]**2)  # (nb_edges_waterline,)

    outputs_per_freq = []
    for i_freq in range(dataset.sizes[freq_type]):
        # One frequency, all the wave directions at once: the expansions have the wave directions in their orders 1
        # and 2, and the shapes in the comments below are their fields.
        ds_at_freq = dataset.isel({freq_type: i_freq})
        omega = float(ds_at_freq.coords["omega"])
        q = Expansion.first_order(motion_amplitude[i_freq, :, :])  # order 1, (nb_dofs,)
        p1_hull, grad_phi1_hull, p1_waterline = _first_order_fields_from_dataset(ds_at_freq, body, mesh, solver, q.coefficients[1])
        grad_p1_hull = 1j * omega * rho * grad_phi1_hull  # order 1, (nb_faces, 3)

        # Hull: Bernoulli + displacement of the pressure to X+u
        u, _ = body._displacement_expansions(q)  # order 1+2, (nb_faces, 3)
        p_expansion = p0 + p1_hull - rho/2 * (grad_phi1_hull * grad_phi1_hull).sum(-1)  # order 0+1+2, (nb_faces,)
        gradient_p_expansion = grad_p0 + grad_p1_hull  # order 0+1, (nb_faces, 3)
        p_displaced = p_expansion + (gradient_p_expansion * u).sum(-1)  # order 0+1+2, (nb_faces,)

        # Waterline: strip between the mean waterline and the free surface, (nb_edges_waterline,)
        eta = p1_waterline / (rho * g)  # free surface elevation, order 1
        u3 = (q.expand_dims(-1) * vertical_dofs_motion_waterline).sum(axis=-2)  # vertical motion of the hull, order 1
        # The second-order motion at the waterline would only multiply order-0 quantities, which are zero there.
        strip_pressure = p1_waterline - rho * g * (eta + u3) / 2 # mean over the height of the strip, order 1
        strip_height = eta - u3  # order 1
        strip_width = strip_height / inclination_waterline  # along the (possibly inclined) hull, order 1
        waterline_force_density = strip_pressure * strip_width  # order 2

        force = (body.integrate_pressure_on_moving_body(p_displaced, q)
                 + body.integrate_pressure_on_waterline(waterline_force_density))
        outputs_per_freq.append((force.coefficients[2], p_displaced.coefficients[2],
                                 strip_pressure.coefficients[1], strip_height.coefficients[1]))

    return _near_field_output_dataset(rao, freq_type, body, mesh, outputs_per_freq, output_pressure)


def _first_order_fields_from_dataset(dataset, body, mesh, solver, q1):
    """First-order fields of the freely floating body at one frequency, for each wave direction:
    incident, diffracted, and radiated with the motion `q1[i_dir, i_dof]`.
    `mesh` is the hull mesh of `body`, merged (for the waterline of symmetric meshes).

    Parameters
    ----------
    q1: array of shape (nb_dir, nb_dofs)
        Complex amplitudes of the dofs for each wave direction.

    Returns
    -------
    First-order DifferenceFrequencyExpansion (order 1 of shape (nb_dir, *field)):
    p1_hull, field (nb_faces,)
        Pressure at the faces centers of the hull.
    grad_phi1_hull, field (nb_faces, 3)
        Gradient of the potential at the faces centers of the hull.
    p1_waterline, field (nb_edges_waterline,)
        Pressure at the middles of the edges of the waterline.
    """
    dofs = list(body.dofs)
    omega = float(dataset.coords["omega"])
    rho = float(dataset.coords["rho"])
    g = float(dataset.coords["g"])
    water_depth = float(dataset.coords["water_depth"])
    p1_hull = (dataset["pressure_hull_incident"] + dataset["pressure_hull_diffraction"]).transpose("wave_direction", "hull_face").values
    p1_hull = p1_hull + q1 @ dataset["pressure_hull_radiation"].sel(radiating_dof=dofs).transpose("radiating_dof", "hull_face").values

    # Diffracted and radiated fields, from the total source distribution for each wave direction
    # (the fields are linear in the sources).
    sources = (_sources_on_mesh_including_lid(dataset, "diffraction", "wave_direction", body)
               + q1 @ _sources_on_mesh_including_lid(dataset.sel(radiating_dof=dofs), "radiation", "radiating_dof", body))
    if not np.all(np.isfinite(sources)):
        raise ValueError(f"The dataset contains missing (NaN) source distributions at omega={omega}, needed for the near field mean drift force.")
    gf_params = dict(free_surface=0.0, water_depth=water_depth, wavenumber=float(dataset.coords["wavenumber"]))
    grad_phi1_hull = solver._compute_potential_gradient_from_sources(mesh, body.mesh_including_lid, sources, **gf_params)
    phi1_waterline = solver._compute_potential_from_sources(mesh.waterline_edges_centers, body.mesh_including_lid, sources, **gf_params)
    p1_waterline = 1j * omega * rho * phi1_waterline

    # Incident waves
    for i_dir, wave_direction in enumerate(dataset.coords["wave_direction"].values):
        problem = DiffractionProblem(body=body, wave_direction=wave_direction, omega=omega, rho=rho, g=g, water_depth=water_depth)
        grad_phi1_hull[i_dir] += airy_waves_velocity(mesh, problem)
        p1_waterline[i_dir] += airy_waves_pressure(mesh.waterline_edges_centers, problem)

    return Expansion.first_order(p1_hull), Expansion.first_order(grad_phi1_hull), Expansion.first_order(p1_waterline)


def _check_inputs_for_near_field(dataset, body, solver):
    """Shallow consistency checks of the dataset input of `compute_near_field_mean_drift_force`."""
    if dataset.attrs.get("method") == "direct":
        raise ValueError("The near field mean drift force needs the source distributions of the indirect method, "
                         "but the dataset has been computed with the direct method.")
    missing = [name for name in ["sources_diffraction", "sources_radiation"] if name not in dataset]
    if body.lid_mesh is not None:
        missing += [name for name in ["sources_lid_diffraction", "sources_lid_radiation"] if name not in dataset]
    if len(missing) > 0:
        raise ValueError(f"The dataset does not contain {missing}, needed for the near field mean drift force.\n"
                         "Compute the dataset with the indirect method and `keep_details=True`, "
                         "e.g. `solver.fill_dataset(test_matrix, body, keep_details=True)`.")
    for name in ["rho", "g", "water_depth", "forward_speed"]:
        if dataset.coords[name].ndim != 0:
            raise NotImplementedError(f"The near field mean drift force only supports datasets with a single value of {name}.")
    if float(dataset.coords["forward_speed"]) != 0.0:
        raise NotImplementedError("The near field mean drift force is not implemented with forward speed.")
    for name in ["omega", "wavenumber"]:
        if name not in dataset.coords:
            raise ValueError(f"The near field mean drift force needs the coordinate {name} in the dataset.")
    if body.lid_mesh is None and "lid_face" in dataset.dims:
        raise ValueError(f"The dataset has been computed with a lid, but {body} has no lid mesh.")
    if not isinstance(body, (FloatingBody, Multibody)):
        raise NotImplementedError(
            f"The near field mean drift force is only implemented for a FloatingBody or a Multibody, not for {body}."
        )
    if not all(isinstance(dof, AbstractDof) for dof in body.dofs.values()):
        raise NotImplementedError(
            "The near field mean drift force is only implemented for dofs defined as AbstractDof "
            "(e.g. TranslationDof, RotationDof, CustomDof), not as arrays."
        )
    if dataset.sizes["hull_face"] != body.mesh.nb_faces or (body.lid_mesh is not None and dataset.sizes["lid_face"] != body.lid_mesh.nb_faces):
        raise ValueError(f"The faces of the dataset do not match the mesh of {body}.")
    # With `mesh=True`, the dataset also stores the faces centers.
    if "mesh_hull_faces_centers" in dataset and not np.allclose(dataset["mesh_hull_faces_centers"].values, body.mesh.faces_centers):
        raise ValueError(f"The faces of the dataset do not match the mesh of {body}.")
    if "mesh_lid_faces_centers" in dataset and not np.allclose(dataset["mesh_lid_faces_centers"].values, body.lid_mesh.faces_centers):
        raise ValueError(f"The faces of the lid of the dataset do not match the lid mesh of {body}.")
    # Only the Green function matters to evaluate the potential from the source distributions.
    green_function = getattr(solver.engine, "green_function", None)
    green_function_settings = green_function.exportable_settings if green_function is not None else {}
    for key, value in green_function_settings.items():
        if key in dataset.attrs and str(dataset.attrs[key]) != str(value):
            LOG.warning(f"The setting {key} of the solver ({value}) differs from the one used to compute the dataset "
                        f"({dataset.attrs[key]}): the potential evaluated from the source distributions might be inconsistent.")


def _align_coordinates(motion_amplitude, dataset, body):
    """Restrict the inputs of `compute_near_field_mean_drift_force` to the frequencies and wave directions
    common to `motion_amplitude` and `dataset`, and to the dofs common to them and `body`.
    A warning is displayed when some values are dropped, and an error when nothing is left.

    Returns
    -------
    motion_amplitude: DataArray of dims (freq_type, "wave_direction", "radiating_dof")
    dataset: Dataset, at the same frequencies and wave directions, in the same order
    body: FloatingBody, with only the common dofs
    freq_type: str
    """
    frequency_keys = set(motion_amplitude.dims) & set(dataset.dims) & {'omega', 'freq', 'period', 'wavelength', 'wavenumber'}
    if len(frequency_keys) == 1:
        freq_type = list(frequency_keys)[0]
    else:
        raise ValueError("Cannot identify dimension holding the frequency in the inputs of `compute_near_field_mean_drift_force`")

    def intersection(name, values, *others):
        # Keep the order of `values`
        common = [v for v in values if all(v in other for other in others)]
        if len(common) == 0:
            raise ValueError(f"The inputs of `compute_near_field_mean_drift_force` have no common {name}.")
        if any(len(common) != len(set(other)) for other in (values, *others)):
            LOG.warning(f"The inputs of `compute_near_field_mean_drift_force` do not have the same {name}: "
                        f"the near field mean drift force is computed only for the common ones.")
        return common

    freqs = intersection(
        freq_type,
        motion_amplitude.coords[freq_type].values,
        dataset.coords[freq_type].values
    )
    wave_dirs = intersection(
        "wave_direction",
        motion_amplitude.coords["wave_direction"].values,
        dataset.coords["wave_direction"].values
    )
    dofs = intersection(
        "radiating_dof",
        list(body.dofs),
        motion_amplitude.coords["radiating_dof"].values,
        dataset.coords["radiating_dof"].values
    )

    all_coords = {freq_type: freqs, "wave_direction": wave_dirs, "radiating_dof": dofs}
    motion_amplitude = motion_amplitude.sel(**all_coords).transpose(freq_type, "wave_direction", "radiating_dof")
    dataset = dataset.sel(**all_coords)
    if len(dofs) < len(body.dofs):
        body = body.with_only_dofs(dofs)
    return motion_amplitude, dataset, body, freq_type


def _sources_on_mesh_including_lid(dataset, kind, leading_dim, body):
    """Source distributions `sources_{kind}` of the dataset (and `sources_lid_{kind}`), on all the faces of
    `body.mesh_including_lid`, as an array of shape (nb_{leading_dim}, nb_faces_including_lid)."""
    hull = dataset[f"sources_{kind}"].transpose(leading_dim, "hull_face").values
    if body.lid_mesh is None:
        return hull
    sources = np.empty((hull.shape[0], body.mesh_including_lid.nb_faces), dtype=complex)
    sources[:, body.hull_mask] = hull
    sources[:, ~body.hull_mask] = dataset[f"sources_lid_{kind}"].transpose(leading_dim, "lid_face").values
    return sources


def _near_field_output_dataset(rao, freq_type, body, mesh, outputs_per_freq, output_pressure):
    """Output dataset of `compute_near_field_mean_drift_force`, from the tuples (force, pressure on the hull,
    pressure and height of the strip at the waterline) of each frequency."""
    force, pressure_hull_mean_drift, pressure_waterline, relative_elevation_waterline = (
            np.stack(output) for output in zip(*outputs_per_freq))
    coords = {
        freq_type: rao.coords[freq_type].values,
        "wave_direction": rao.coords["wave_direction"].values,
        "other_wave_direction": rao.coords["wave_direction"].values,
        }
    output = xr.Dataset(
            {"mean_drift_force_near_field": xr.DataArray(
                data=force,
                dims=[freq_type, "wave_direction", "other_wave_direction", "influenced_dof"],
                coords={**coords, "influenced_dof": list(body.dofs.keys())},
                )
             })

    if output_pressure:
        output["pressure_hull_mean_drift"] = xr.DataArray(
                data=pressure_hull_mean_drift,
                dims=[freq_type, "wave_direction", "other_wave_direction", "hull_face"],
                coords=coords,
                 )
        waterline_coords = {
            freq_type: rao.coords[freq_type].values,
            "wave_direction": rao.coords["wave_direction"].values,
            }
        output["pressure_waterline_total"] = xr.DataArray(
                pressure_waterline,
                dims=[freq_type, "wave_direction", "waterline_edge"],
                coords=waterline_coords,
                 )
        output["relative_elevation_waterline"] = xr.DataArray(
                relative_elevation_waterline,
                dims=[freq_type, "wave_direction", "waterline_edge"],
                coords=waterline_coords,
                 )
        output["mesh_waterline_vertices"] = xr.DataArray(
                mesh.vertices[mesh.edges_waterline, :],
                dims=["waterline_edge", "vertices_of_edge", "space_coordinate"],
                )

    return output
