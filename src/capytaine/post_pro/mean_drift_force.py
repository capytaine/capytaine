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

import numpy as np
import xarray as xr

from capytaine.bem.airy_waves import (
    airy_waves_pressure,
    airy_waves_velocity,
)
from capytaine.bodies.bodies import FloatingBody
from capytaine.bodies.multibodies import Multibody
from capytaine.bodies.dofs import AbstractDof
from capytaine.bem.problems_and_results import DiffractionResult, RadiationResult
from capytaine.tools.perturbation_expansion import DifferenceFrequencyExpansion as Expansion


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

def near_field_mean_drift_force(first_order_motion_amplitude, results, solver, *, output_pressure=False):
    """Compute the mean drift force using the near field formulation.

    The force is the second-order part of the integral of the pressure on the moving body and along
    the waterline, for any dofs (rigid body motions or generalized modes with a gradient of motion).
    All the wave components have the same frequency, so all the terms of the result (for the same or for
    different wave directions) are at difference frequency zero, where the second-order potential does
    not contribute to the pressure.

    Current assumptions: all results are for the same :class:`~capytaine.bodies.bodies.FloatingBody`
    or :class:`~capytaine.bodies.multibodies.Multibody`, the same rho and g.
    Several wavelengths and dofs are taken into account at once.

    Parameters
    ----------
    first_order_motion_amplitude : xarray DataArray
        The motion RAO, as returned by :func:`capytaine.post_pro.rao.rao`.
    results : list of LinearPotentialFlowResult
        The diffraction and radiation results used to compute the RAO, for all the wavelengths
        and wave directions of interest.
    solver : BEMSolver
        The solver used to compute the potential and its gradient on the hull at the requested
        points.
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
    body = results[0].body
    if isinstance(body, Multibody):
        body = body.as_FloatingBody  # Same mesh and dofs (`DofOnSubmesh`), with the methods of FloatingBody.
    if not isinstance(body, FloatingBody):
        raise NotImplementedError(f"The near field mean drift force is only implemented for a FloatingBody or a Multibody, not for {body}.")
    if not all(isinstance(dof, AbstractDof) for dof in body.dofs.values()):
        raise NotImplementedError("The near field mean drift force is only implemented for dofs defined as AbstractDof (e.g. TranslationDof, RotationDof, CustomDof), not as arrays.")
    mesh = body.mesh.merged()
    rho = results[0].rho
    g = results[0].g
    freq_type = results[0].provided_freq_type
    first_order_motion_amplitude = first_order_motion_amplitude.transpose(freq_type, "wave_direction", "radiating_dof")
    omegas = first_order_motion_amplitude.coords["omega"].values
    motion_amplitude = first_order_motion_amplitude.sel(radiating_dof=list(body.dofs)).values # motion_amplitude[i_freq, i_dir, i_dof]

    # First-order amplitudes
    pressure_hull = total_pressure_on_hull(mesh, results, first_order_motion_amplitude).values # pressure_hull[i_freq, i_dir, i_face]
    gradient_potential = total_potential_gradient(solver, mesh, results, first_order_motion_amplitude).values # gradient_potential[i_freq, i_dir, i_face, xyz]
    hydrodynamic_pressure_waterline = total_pressure(solver, mesh.middle_edges_waterline, results, first_order_motion_amplitude).values # hydrodynamic_pressure_waterline[i_freq, i_dir, i_vertex_waterline]
    dofs_motion_waterline = body._dofs_on_waterline() # dofs_motion_waterline[i_dof, i_vertex_waterline, xyz]
    vertical_motion_waterline = np.sum(motion_amplitude[..., None] * dofs_motion_waterline[None, None, :, :, 2], axis=2) # vertical_motion_waterline[i_freq, i_dir, i_vertex_waterline]
    free_surface_elevation_waterline = hydrodynamic_pressure_waterline / (rho * g)
    # Pressure averaged over the height of the strip between the mean waterline and the free surface,
    # height of this strip, and its width along the (possibly inclined) hull.
    pressure_waterline = hydrodynamic_pressure_waterline - rho * g * (free_surface_elevation_waterline + vertical_motion_waterline)/2
    relative_elevation_waterline = free_surface_elevation_waterline - vertical_motion_waterline
    strip_width_waterline = relative_elevation_waterline / np.sqrt(1 - mesh.normals_edges_waterline[:, 2]**2)

    # Hydrostatic pressure on the hull (order 0)
    hydrostatic_pressure = - rho * g * mesh.faces_centers[:, 2] # hydrostatic_pressure[i_face]
    gradient_hydrostatic_pressure = np.array([0.0, 0.0, - rho * g])

    nb_freq, nb_dir = motion_amplitude.shape[:2]
    F = np.empty((nb_freq, nb_dir, nb_dir, len(body.dofs)), dtype=complex)
    p2 = np.empty((nb_freq, nb_dir, nb_dir, mesh.nb_faces), dtype=complex)
    for w in range(nb_freq):
        u, _ = body._displacement_expansions(motion_amplitude[w]) # field (nb_faces, 3)
        velocity = Expansion.first_order(gradient_potential[w]) # field (nb_faces, 3)
        gradient_pressure = gradient_hydrostatic_pressure + Expansion.first_order(1j * omegas[w] * rho * gradient_potential[w]) # field (nb_faces, 3)
        # Bernoulli (without second-order potential) at the displaced point X + u
        pressure = hydrostatic_pressure + Expansion.first_order(pressure_hull[w]) - rho/2 * (velocity * velocity).sum(-1) + (gradient_pressure * u).sum(-1) # field (nb_faces,)
        # The second-order motion at the waterline would only multiply order-0 quantities, which are zero at the waterline.
        waterline_term = Expansion.first_order(pressure_waterline[w]) * Expansion.first_order(strip_width_waterline[w]) # field (nb_waterline_segments,)
        force = body.integrate_pressure_on_moving_body(pressure, motion_amplitude[w]) + body.integrate_pressure_on_waterline(waterline_term)
        F[w] = force.coefficients[2]
        p2[w] = pressure.coefficients[2]

    coords = {
        freq_type: first_order_motion_amplitude.coords[freq_type].values,
        "wave_direction_k": first_order_motion_amplitude.coords["wave_direction"].values,
        "wave_direction_l": first_order_motion_amplitude.coords["wave_direction"].values,
        }
    dataset = xr.Dataset(
            {"near_field_mean_drift_force": xr.DataArray(
                data=F,
                dims=[freq_type, "wave_direction_k", "wave_direction_l", "influenced_dof"],
                coords={**coords, "influenced_dof": list(body.dofs.keys())},
                )
             })

    if output_pressure:
        dataset["pressure_hull_mean_drift"] = xr.DataArray(
                data=p2,
                dims=[freq_type, "wave_direction_k", "wave_direction_l", "hull_face"],
                coords=coords,
                 )
        waterline_coords = {
            freq_type: first_order_motion_amplitude.coords[freq_type].values,
            "wave_direction": first_order_motion_amplitude.coords["wave_direction"].values,
            }
        dataset["pressure_waterline_total"] = xr.DataArray(
                pressure_waterline,
                dims=[freq_type, "wave_direction", "waterline_segment"],
                coords=waterline_coords,
                 )
        dataset["relative_elevation_waterline"] = xr.DataArray(
                relative_elevation_waterline,
                dims=[freq_type, "wave_direction", "waterline_segment"],
                coords=waterline_coords,
                 )
        dataset["mesh_waterline"] = xr.DataArray(
                mesh.vertices[mesh.edges_waterline, :],
                dims=["waterline_segment", "point_per_segment", "space_coordinate"],
                )

    return dataset

def total_potential_gradient(solver, mesh, results, rao):
    freq_type = results[0].provided_freq_type
    potential_gradient = xr.DataArray(
        data=np.zeros((rao.sizes[freq_type], rao.sizes["wave_direction"], mesh.nb_faces, 3), dtype=complex),
        coords={
            freq_type: rao.coords[freq_type],
            "wave_direction": rao.coords["wave_direction"],
            "mesh_face": np.arange(mesh.nb_faces),
            "xyz": np.arange(3),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            potential_gradient.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += airy_waves_velocity(mesh, res) # incident
            potential_gradient.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += solver._compute_potential_gradient(mesh, res) # diffracted

        elif isinstance(res, RadiationResult):
            rao_rad = rao.sel({freq_type: getattr(res, freq_type), "radiating_dof": res.radiating_dof}).values
            potential_gradient.loc[{freq_type: getattr(res, freq_type)}] += rao_rad[:, None, None] * solver._compute_potential_gradient(mesh, res)[None, :, :] # radiated
            # potential_gradient[i_freq, i_dir, i_face, xyz] = rao_rad[i_dir] * _compute_potential_gradient[i_face, xyz]

    return potential_gradient

def total_pressure(solver, vertices_middle, results, rao):
    nb_vertices_waterline = np.shape(vertices_middle)[0]
    freq_type = results[0].provided_freq_type
    pressure = xr.DataArray(
        data=np.zeros((rao.sizes[freq_type], rao.sizes["wave_direction"], nb_vertices_waterline), dtype=complex),
        coords={
            freq_type: rao.coords[freq_type],
            "wave_direction": rao.coords["wave_direction"],
            "waterline_segment": np.arange(nb_vertices_waterline),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            pressure.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += airy_waves_pressure(vertices_middle, res) # incident
            pressure.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += solver.compute_pressure(vertices_middle, res) # diffracted

        elif isinstance(res, RadiationResult):
            rao_rad = rao.sel({freq_type: getattr(res, freq_type), "radiating_dof": res.radiating_dof}).values
            pressure.loc[{freq_type: getattr(res, freq_type)}] += rao_rad[:, None] * solver.compute_pressure(vertices_middle, res)[None, :] # radiated
            # pressure[i_freq, i_dir, i_edge_waterline] = rao_rad[i_dir] * compute_pressure[i_edge_waterline]

    return pressure

def total_pressure_on_hull(mesh, results, rao):
    freq_type = results[0].provided_freq_type
    pressure = xr.DataArray(
        data=np.zeros((rao.sizes[freq_type], rao.sizes["wave_direction"], mesh.nb_faces), dtype=complex),
        coords={
            freq_type: rao.coords[freq_type],
            "wave_direction": rao.coords["wave_direction"],
            "hull_face": np.arange(mesh.nb_faces),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            pressure.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += airy_waves_pressure(mesh, res) # incident
            pressure.loc[{freq_type: getattr(res, freq_type), "wave_direction": res.wave_direction}] += res.pressure_on_hull # diffracted

        elif isinstance(res, RadiationResult):
            rao_rad = rao.sel({freq_type: getattr(res, freq_type), "radiating_dof": res.radiating_dof}).values
            pressure.loc[{freq_type: getattr(res, freq_type)}] += rao_rad[:, None] * res.pressure_on_hull[None, :] # radiated
            # pressure[i_freq, i_dir, i_face] = rao_rad[i_dir] * pressure_on_hull[i_face]

    return pressure
