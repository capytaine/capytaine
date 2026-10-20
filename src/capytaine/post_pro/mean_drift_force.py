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

def far_field_mean_drift_force(X, dataset):
    """Compute the mean drift forces using far field formulation.
    Note that the forces are proportional to the square of the wave amplitude, but this implementation
    does not take into account the wave amplitude.

    Parameters
    ----------
    X : xarray DataArray
        The motion RAO.
    dataset : xarray Dataset
        This function supposes that variables named 'kochin_diffraction' and 'kochin_radiation' are in the dataset.

    Returns
    -------
    xarray Dataset
        The horizontal mean drift forces, depending on omega and the wave direction.
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

    H_rad_tot = sum(H_rad.sel(radiating_dof=d)*X.sel(radiating_dof=d) for d in X.radiating_dof)
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

    dims = [X.dims[0], "wave_direction_k", "wave_direction_l"]
    coords = {
        X.dims[0]: X.coords[X.dims[0]].values,
        "wave_direction_k": X.coords["wave_direction"].values,
        "wave_direction_l": X.coords["wave_direction"].values,
    }

    freq = X.sizes[X.dims[0]]
    nb_dir = X.sizes["wave_direction"]
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

    Fx = xr.DataArray(
        data=data_x,
        dims=dims,
        coords=coords,
        name='drift_force_surge'
        )

    Fy = xr.DataArray(
        data=data_y,
        coords=coords,
        dims=dims,
        name='drift_force_sway'
        )

    Mz = xr.DataArray(
        data=data_z,
        coords=coords,
        dims=dims,
        name='drift_force_yaw'
        )
    return xr.Dataset({Fx.name: Fx, Fy.name: Fy, Mz.name: Mz})


def _merge_far_field_mean_drift_variables(dataset):
    """Merge the three drift force components into a single variable with an influenced_dof dimension.

    Parameters
    ----------
    dataset : xarray Dataset
        The dataset returned by far_field_mean_drift_force, containing 'drift_force_surge',
        'drift_force_sway', and 'drift_force_yaw' variables.

    Returns
    -------
    xarray Dataset
        A dataset with a single 'drift_force' variable having an additional 'influenced_dof' dimension
        with coordinates ["Surge", "Sway", "Yaw"].
    """
    # Stack the three force components along a new dimension
    drift_force = xr.concat(
        [dataset['drift_force_surge'], dataset['drift_force_sway'], dataset['drift_force_yaw']],
        dim=xr.DataArray(['Surge', 'Sway', 'Yaw'], dims=['influenced_dof'], name='influenced_dof')
    )
    drift_force.name = 'drift_force'

    return xr.Dataset({'drift_force': drift_force})


####################################################################################################
# Near field
####################################################################################################

def near_field_mean_drift_force(first_order_motion_amplitude, dataset, body, *, solver=None, output_pressure=False):
    """Compute the mean drift force using the near field formulation.

    The force is the second-order part of the integral of the pressure on the moving body and along
    the waterline, for any dofs (rigid body motions or generalized modes with a gradient of motion).
    All the wave components have the same frequency, so all the terms of the result (for the same or for
    different wave directions) are at difference frequency zero, where the second-order potential does
    not contribute to the pressure.

    The first-order fields are evaluated from a dataset computed with the indirect method and
    ``keep_details=True`` (e.g. ``solver.fill_dataset(test_matrix, body, keep_details=True)``),
    which stores the pressures and the source distributions on the body.
    Several wavelengths and dofs are taken into account at once.
    Only datasets with a single value of rho, g and water depth, and without forward speed, are supported.

    Parameters
    ----------
    first_order_motion_amplitude : xarray DataArray
        The motion RAO, as returned by :func:`capytaine.post_pro.rao.rao`.
    dataset : xarray Dataset
        The dataset of the diffraction and radiation problems used to compute the RAO, for all the
        wavelengths and wave directions of interest, with the source distributions on the body
        (``sources_diffraction``, ``sources_radiation``, and on the lid if any).
    body : FloatingBody or Multibody
        The body of the dataset.
    solver : BEMSolver, optional
        The solver evaluating the potential and its gradient from the source distributions.
        It should use the same Green function as the one used to compute the dataset (a warning
        is raised when the settings stored in the dataset differ).
        Default: ``BEMSolver()``.
    output_pressure : bool, optional
        If True, also return the second-order pressure field on the hull and the first-order
        fields on the waterline from which the waterline term of the force is computed.
        Default: False.

    Returns
    -------
    xarray Dataset
        A dataset with a ``near_field_mean_drift_force`` variable of shape
        ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_dofs)``.
        If ``output_pressure`` is True, the dataset also contains the following variables:

        - ``pressure_hull_mean_drift``, of shape
          ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_hull_faces)``: the mean drift
          pressure field on the hull. It has the same normalization as
          ``near_field_mean_drift_force``: its diagonal is the mean pressure in a single wave of
          unit amplitude. It does not include the waterline term.
          For a fixed body, the hull term of the force is the integral of this pressure on the mean hull.
          For a moving body, it also includes the first-order pressure times the first-order
          variations of the normals and of the motion of the dofs on the moving hull.
        - ``pressure_waterline_total``, of shape
          ``(nb_freq, nb_wave_direction, nb_waterline_segments)``: the complex first-order
          pressure ``p`` in the strip between the mean waterline and the free surface, at the
          middle of each segment of the waterline, for a wave of unit amplitude.
          It is the total pressure on the freely floating body: incident, diffracted and
          radiated (with the motion given by ``first_order_motion_amplitude``) waves, as well as the hydrostatic pressure
          averaged over the height of the strip.
        - ``relative_elevation_waterline``, of the same shape: the complex first-order relative
          elevation ``h`` of the free surface with respect to the body, at the same points, for
          a wave of unit amplitude.
        - ``mesh_waterline``, of shape ``(nb_waterline_segments, 2, 3)``: the coordinates of
          the two ends of each segment of the waterline.

        For the waves ``k`` and ``l``, the waterline term of ``near_field_mean_drift_force`` is
        the integral along the waterline of ``(p_k conj(w_l) + conj(p_l) w_k)/4 * (-n·δu_j)`` for the dof ``j``,
        where ``w = h / sqrt(1 - n_z²)`` is the width of the strip along the hull and ``n`` the normal of the
        face containing the segment of the waterline.
    """
    solver, floating_body, freq_type, motion_amplitude, dataset = _prepare_near_field_inputs(
            first_order_motion_amplitude, dataset, body, solver)
    mesh = floating_body.mesh.merged()  # The waterline data of symmetric meshes is computed on the merged mesh.
    rho = float(dataset.coords["rho"])
    g = float(dataset.coords["g"])

    # Still water (order 0): hydrostatic pressure on the hull and geometry of the waterline
    p0 = -rho * g * mesh.faces_centers[:, 2]  # (nb_faces,)
    grad_p0 = np.array([0.0, 0.0, -rho * g])
    vertical_dofs_motion_waterline = floating_body._dofs_on_waterline()[:, :, 2]  # (nb_dofs, nb_edges_waterline)
    inclination_waterline = np.sqrt(1 - mesh.normals_edges_waterline[:, 2]**2)  # (nb_edges_waterline,)

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
        u, _ = floating_body._displacement_expansions(q)  # order 1+2, (nb_faces, 3)
        p_expansion = p0 + p1_hull - rho/2 * (grad_phi1_hull * grad_phi1_hull).sum(-1)  # order 0+1+2, (nb_faces,)
        gradient_p_expansion = grad_p0 + grad_p1_hull  # order 0+1, (nb_faces, 3)
        p_displaced = p_expansion + (gradient_p_expansion * u).sum(-1)  # order 0+1+2, (nb_faces,)

        # Waterline: strip between the mean waterline and the free surface, (nb_edges_waterline,)
        # The second-order motion at the waterline would only multiply order-0 quantities, which are zero there.
        eta = p1_waterline / (rho * g)  # free surface elevation, order 1
        u3 = (q.expand_dims(-1) * vertical_dofs_motion_waterline).sum(axis=-2)  # vertical motion of the hull, order 1
        strip_pressure = p1_waterline - rho * g * (eta + u3) / 2 # mean over the height of the strip, order 1
        strip_height = eta - u3  # order 1
        strip_width = strip_height / inclination_waterline  # along the (possibly inclined) hull, order 1
        waterline_force_density = strip_pressure * strip_width  # order 2

        force = (floating_body.integrate_pressure_on_moving_body(p_displaced, q)
                 + floating_body.integrate_pressure_on_waterline(waterline_force_density))
        outputs_per_freq.append((force.coefficients[2], p_displaced.coefficients[2],
                                 strip_pressure.coefficients[1], strip_height.coefficients[1]))

    return _near_field_output_dataset(first_order_motion_amplitude, freq_type, body, mesh, outputs_per_freq, output_pressure)


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
    p1_hull = (dataset["pressure_Froude_Krylov"] + dataset["pressure_diffraction"]).transpose("wave_direction", "hull_face").values
    p1_hull = p1_hull + q1 @ dataset["pressure_radiation"].sel(radiating_dof=dofs).transpose("radiating_dof", "hull_face").values

    # Diffracted and radiated fields, from the total source distribution for each wave direction
    # (the fields are linear in the sources).
    sources = (_sources_on_mesh_including_lid(dataset, "diffraction", "wave_direction", body)
               + q1 @ _sources_on_mesh_including_lid(dataset.sel(radiating_dof=dofs), "radiation", "radiating_dof", body))
    if not np.all(np.isfinite(sources)):
        raise ValueError(f"The dataset contains missing (NaN) source distributions at omega={omega}, needed for the near field mean drift force.")
    gf_params = dict(free_surface=0.0, water_depth=water_depth, wavenumber=float(dataset.coords["wavenumber"]))
    grad_phi1_hull = solver._compute_potential_gradient_from_sources(mesh, body.mesh_including_lid, sources, **gf_params)
    phi1_waterline = solver._compute_potential_from_sources(mesh.middle_edges_waterline, body.mesh_including_lid, sources, **gf_params)
    p1_waterline = 1j * omega * rho * phi1_waterline

    # Incident waves
    for i_dir, wave_direction in enumerate(dataset.coords["wave_direction"].values):
        problem = DiffractionProblem(body=body, wave_direction=wave_direction, omega=omega, rho=rho, g=g, water_depth=water_depth)
        grad_phi1_hull[i_dir] += airy_waves_velocity(mesh, problem)
        p1_waterline[i_dir] += airy_waves_pressure(mesh.middle_edges_waterline, problem)

    return Expansion.first_order(p1_hull), Expansion.first_order(grad_phi1_hull), Expansion.first_order(p1_waterline)


def _prepare_near_field_inputs(first_order_motion_amplitude, dataset, body, solver):
    """Checks and defaults of the inputs of `near_field_mean_drift_force`.

    Returns
    -------
    solver: BEMSolver
    floating_body: FloatingBody
        `body`, or `body.as_FloatingBody` for a Multibody (same mesh and dofs, `DofOnSubmesh` for a Multibody),
        whose methods compute the second-order part. The first-order fields are evaluated with `body` itself,
        whose mesh including lid and hull mask are those of the dataset.
    freq_type: str
        The frequency dimension of `first_order_motion_amplitude`.
    motion_amplitude: array of shape (nb_freq, nb_dir, nb_dofs)
    dataset: xarray Dataset
        The dataset at the frequencies and wave directions of `first_order_motion_amplitude`, in the same order.
    """
    if solver is None:
        from capytaine.bem.solver import BEMSolver  # Not at the top of the module, to avoid a circular import.
        solver = BEMSolver()
    floating_body = body.as_FloatingBody if isinstance(body, Multibody) else body
    if not isinstance(floating_body, FloatingBody):
        raise NotImplementedError(f"The near field mean drift force is only implemented for a FloatingBody or a Multibody, not for {body}.")
    if not all(isinstance(dof, AbstractDof) for dof in body.dofs.values()):
        raise NotImplementedError("The near field mean drift force is only implemented for dofs defined as AbstractDof (e.g. TranslationDof, RotationDof, CustomDof), not as arrays.")
    _check_dataset_for_near_field(dataset, body, solver)
    freq_type = _check_motion_amplitude_for_near_field(first_order_motion_amplitude, dataset, body)

    first_order_motion_amplitude = first_order_motion_amplitude.transpose(freq_type, "wave_direction", "radiating_dof")
    motion_amplitude = first_order_motion_amplitude.sel(radiating_dof=list(body.dofs)).values
    dataset = dataset.sel({freq_type: first_order_motion_amplitude.coords[freq_type].values,
                           "wave_direction": first_order_motion_amplitude.coords["wave_direction"].values})
    return solver, floating_body, freq_type, motion_amplitude, dataset


def _check_dataset_for_near_field(dataset, body, solver):
    """Shallow consistency checks of the dataset input of `near_field_mean_drift_force`."""
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
    if dataset.sizes["hull_face"] != body.mesh.nb_faces or (body.lid_mesh is not None and dataset.sizes["lid_face"] != body.lid_mesh.nb_faces):
        raise ValueError(f"The faces of the dataset do not match the mesh of {body}.")
    # With `mesh=True`, the dataset also stores the faces centers.
    if "mesh_faces_center" in dataset and not np.allclose(dataset["mesh_faces_center"].values, body.mesh.faces_centers):
        raise ValueError(f"The faces of the dataset do not match the mesh of {body}.")
    if "lid_mesh_faces_center" in dataset and not np.allclose(dataset["lid_mesh_faces_center"].values, body.lid_mesh.faces_centers):
        raise ValueError(f"The faces of the lid of the dataset do not match the lid mesh of {body}.")
    missing_dofs = set(body.dofs) - set(dataset.coords["radiating_dof"].values)
    if len(missing_dofs) > 0:
        raise ValueError(f"The dataset does not contain the radiation problems of the dofs {sorted(missing_dofs)}.")
    # Only the Green function matters to evaluate the potential from the source distributions.
    green_function = getattr(solver.engine, "green_function", None)
    green_function_settings = green_function.exportable_settings if green_function is not None else {}
    for key, value in green_function_settings.items():
        if key in dataset.attrs and str(dataset.attrs[key]) != str(value):
            LOG.warning(f"The setting {key} of the solver ({value}) differs from the one used to compute the dataset "
                        f"({dataset.attrs[key]}): the potential evaluated from the source distributions might be inconsistent.")


def _check_motion_amplitude_for_near_field(first_order_motion_amplitude, dataset, body):
    """Check that `first_order_motion_amplitude` matches the dataset, and return the name of its frequency dimension."""
    freq_types = [d for d in first_order_motion_amplitude.dims if d not in {"wave_direction", "radiating_dof"}]
    if len(freq_types) != 1 or freq_types[0] not in dataset.coords:
        raise ValueError("The first order motion amplitude should have the dimensions 'wave_direction', 'radiating_dof' "
                         f"and a frequency dimension of the dataset, got {first_order_motion_amplitude.dims}.")
    freq_type, = freq_types
    for name in [freq_type, "wave_direction"]:
        missing = set(first_order_motion_amplitude.coords[name].values) - set(dataset.coords[name].values)
        if len(missing) > 0:
            raise ValueError(f"The values {sorted(missing)} of {name} of the first order motion amplitude are not in the dataset.")
    missing_dofs = set(body.dofs) - set(first_order_motion_amplitude.coords["radiating_dof"].values)
    if len(missing_dofs) > 0:
        raise ValueError(f"The first order motion amplitude does not contain the dofs {sorted(missing_dofs)}.")
    return freq_type


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


def _near_field_output_dataset(first_order_motion_amplitude, freq_type, body, mesh, outputs_per_freq, output_pressure):
    """Output dataset of `near_field_mean_drift_force`, from the tuples (force, pressure on the hull,
    pressure and height of the strip at the waterline) of each frequency."""
    force, pressure_hull_mean_drift, pressure_waterline, relative_elevation_waterline = (
            np.stack(output) for output in zip(*outputs_per_freq))
    coords = {
        freq_type: first_order_motion_amplitude.coords[freq_type].values,
        "wave_direction_k": first_order_motion_amplitude.coords["wave_direction"].values,
        "wave_direction_l": first_order_motion_amplitude.coords["wave_direction"].values,
        }
    output = xr.Dataset(
            {"near_field_mean_drift_force": xr.DataArray(
                data=force,
                dims=[freq_type, "wave_direction_k", "wave_direction_l", "influenced_dof"],
                coords={**coords, "influenced_dof": list(body.dofs.keys())},
                )
             })

    if output_pressure:
        output["pressure_hull_mean_drift"] = xr.DataArray(
                data=pressure_hull_mean_drift,
                dims=[freq_type, "wave_direction_k", "wave_direction_l", "hull_face"],
                coords=coords,
                 )
        waterline_coords = {
            freq_type: first_order_motion_amplitude.coords[freq_type].values,
            "wave_direction": first_order_motion_amplitude.coords["wave_direction"].values,
            }
        output["pressure_waterline_total"] = xr.DataArray(
                pressure_waterline,
                dims=[freq_type, "wave_direction", "waterline_segment"],
                coords=waterline_coords,
                 )
        output["relative_elevation_waterline"] = xr.DataArray(
                relative_elevation_waterline,
                dims=[freq_type, "wave_direction", "waterline_segment"],
                coords=waterline_coords,
                 )
        output["mesh_waterline"] = xr.DataArray(
                mesh.vertices[mesh.edges_waterline, :],
                dims=["waterline_segment", "point_per_segment", "space_coordinate"],
                )

    return output
