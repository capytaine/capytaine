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
    airy_waves_free_surface_elevation,
    airy_waves_velocity,
    froude_krylov_force,
)
from capytaine.bodies.dofs import TranslationDof, RotationDof
from capytaine.meshes.geometry import compute_faces_normals
from capytaine.bem.problems_and_results import DiffractionResult, RadiationResult


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

def near_field_mean_drift_force(rao, results, solver):
    """
    Current assumptions: all results are for the same body, the same frequency and the same rho and g.
    Several wavelengths and dofs are taken into account at once.
    """
    body = results[0].body
    mesh = body.mesh
    rho = results[0].rho
    g = results[0].g
    omega = rao.coords["omega"].values
    nb_dir = rao.sizes["wave_direction"]
    zero_block = np.zeros((nb_dir, 3, 3))
    zero_33 = zero_block[0, ...]

    translation = rao.sel(radiating_dof=[dof_name for dof_name, dof in body.dofs.items() if isinstance(dof, TranslationDof)]).values[0, :, :] # translation[i_dir, xyz]
    rotation = rao.sel(radiating_dof=[dof_name for dof_name, dof in body.dofs.items() if isinstance(dof, RotationDof)]).values[0, :, :] # rotation[i_dir, xyz]
    translation_matrix = np.block([[zero_block, zero_block], [skew_matrix(translation), zero_block]])
    rotation_matrix = np.block([[skew_matrix(rotation), zero_block], [zero_block, skew_matrix(rotation)]])

    z0 = np.tile(mesh.faces_centers[:, -1], (nb_dir, 1))
    motion = motion_order1(mesh.faces_centers, body, rao) # motion[i_dir, i_face, xyz]
    z1 = motion[:, :, -1] # z1[i_dir, i_face]
    forces_order0 = hydrostatics_forces(body, rao, rho, g, z0)[0, ...] # same value for all wave directions at order 0
    forces_order1 = hydrodynamics_forces_order1(results, rao) + hydrostatics_forces(body, rao, rho, g, z1) # forces_order1[i_dir, influenced_dof]
    rotation_forces_order0 = np.matvec(rotation_matrix, forces_order0) # rotation_forces_order0[i_dir, influenced_dof]
    all_forces_order1 = forces_order1 + rotation_forces_order0 # all_forces_order1[i_dir, influenced_dof]

    gradient_potential = total_potential_gradient(solver, mesh, results, rao) # gradient_potential[i_dir, i_face, xyz]
    edges_waterline = mesh.edges_waterline
    vertices_middle_waterline = (mesh.vertices[edges_waterline[:, 0], :] + mesh.vertices[edges_waterline[:, 1], :]) / 2
    free_surface_elevation = total_free_surface_elevation(solver, vertices_middle_waterline, results, rao) # free_surface_elevation[i_dir, i_vertex_waterline]
    vertical_position = motion_order1(vertices_middle_waterline, body, rao)[:, :, -1] # vertical_position[i_dir, i_vertex_waterline]

    F = np.full((nb_dir, nb_dir, 6), np.nan + 1j*np.nan, dtype=complex)
    for k in range(nb_dir):
        for l in range(k, nb_dir):
            h = transformation_matrix(rao.isel(wave_direction=k)[0, ...], rao.isel(wave_direction=l)[0, ...]) # h[xyz, xyz]
            H = np.block([[h, zero_33], [zero_33, h]])
            product = (np.sum(motion[k, ...] * np.conjugate(-1j * omega * gradient_potential[l, ...]), axis=1) + np.sum(np.conjugate(motion[l, ...]) * -1j * omega * gradient_potential[k, ...], axis=1)) / 2 # product[i_face]
            gradient_potential_square = np.sum(gradient_potential[k, ...] * np.conjugate(gradient_potential[l, ...]), axis=1) # gradient_potential_square[i_face]
            z_order2 = np.matvec(h, mesh.faces_centers)[..., -1] # z_order2[i_face]
            pressure_field = - (product + gradient_potential_square/2 + g * z_order2) # pressure_field[i_face]
            waterline_field = (1/2) * g * (free_surface_elevation - vertical_position)[k, ...] * np.conjugate(free_surface_elevation - vertical_position)[l, ...] # waterline_field[i_vertex_waterline]

            hydrostatics_order2 = np.matvec(H, forces_order0)
            rotation_forces_order1 = (np.matvec(rotation_matrix[k, ...], np.conjugate(forces_order1[l, ...])) + np.matvec(np.conjugate(rotation_matrix[l, ...]), forces_order1[k, ...])) / 2
            translation_moment = (np.matvec(translation_matrix[k, ...], np.conjugate(all_forces_order1[l, ...])) + np.matvec(np.conjugate(translation_matrix[l, ...]), all_forces_order1[k, ...])) / 2
            pressure_hull = body.integrate_pressure(pressure_field)
            pressure_waterline = integrate_pressure_waterline(body, waterline_field)

            F[k, l, :] = rotation_forces_order1 + hydrostatics_order2 + translation_moment + rho * (np.array(list(pressure_hull.values())) + np.array(list(pressure_waterline.values())))
            F[l, k, :] = np.conjugate(F[k, l, :])

    return xr.DataArray(
        data=F/2,
        dims=["wave_direction_k", "wave_direction_l", "influenced_dof"],
        coords={
            "wave_direction_k": rao.coords["wave_direction"].values,
            "wave_direction_l": rao.coords["wave_direction"].values,
            "influenced_dof": list(body.dofs.keys()),
        },
        name="near_field_mean_drift_force"
    )

def motion_order1(points, body, rao):
    return sum(rao.sel(radiating_dof=dof).values[0,:,None,None] * body.dofs[dof].evaluate_motion_at_points(points) for dof in body.dofs)

def total_potential_gradient(solver, mesh, results, rao):
    potential_gradient = xr.DataArray(
        data=np.zeros((rao.sizes["wave_direction"], mesh.nb_faces, 3), dtype=complex),
        coords={
            "wave_direction": rao.coords["wave_direction"],
            "mesh_face": np.arange(mesh.nb_faces),
            "xyz": np.arange(3),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            potential_gradient.loc[{"wave_direction": res.wave_direction}] += airy_waves_velocity(mesh, res) # incident
            potential_gradient.loc[{"wave_direction": res.wave_direction}] += solver._compute_potential_gradient(mesh, res) # diffracted

        elif isinstance(res, RadiationResult):
            rao_rad = rao.sel({"radiating_dof": res.radiating_dof}).values[0, :]
            potential_gradient += rao_rad[:, None, None] * solver._compute_potential_gradient(mesh, res)[None, :, :] # radiated
            # potential_gradient[i_dir, i_face, xyz] = rao_rad[i_dir] * _compute_potential_gradient[i_face, xyz]

    return potential_gradient

def total_free_surface_elevation(solver, vertices_middle, results, rao):
    nb_vertices_waterline = np.shape(vertices_middle)[0]
    free_surface_elevation = xr.DataArray(
        data=np.zeros((rao.sizes["wave_direction"], nb_vertices_waterline), dtype=complex),
        coords={
            "wave_direction": rao.coords["wave_direction"],
            "vertices_waterline": np.arange(nb_vertices_waterline),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            free_surface_elevation.loc[{"wave_direction": res.wave_direction}] += airy_waves_free_surface_elevation(vertices_middle, res) # incident
            free_surface_elevation.loc[{"wave_direction": res.wave_direction}] += solver.compute_free_surface_elevation(vertices_middle, res) # diffracted

        elif isinstance(res, RadiationResult):
            rao_rad = rao.sel({"radiating_dof": res.radiating_dof}).values[0, :]
            free_surface_elevation += rao_rad[:, None] * solver.compute_free_surface_elevation(vertices_middle, res)[None, :] # radiated
            # free_surface_elevation[i_dir, i_edge_waterline] = rao_rad[i_dir] * compute_free_surface_elevation[i_edge_waterline]

    return free_surface_elevation

def hydrostatics_forces(body, rao, rho, g, z):
    force = np.zeros((rao.sizes["wave_direction"], rao.sizes["radiating_dof"]), dtype=z.dtype)
    for idx, _ in enumerate(rao.coords["wave_direction"].values):
        results = body.integrate_pressure(-g * z[idx, :])
        force[idx, :] = rho * np.array(list(results.values()))
    return force


def hydrodynamics_forces_order1(results, rao):
    force = np.zeros(list(rao.sizes.values()), dtype=complex)

    for res in results:
        if isinstance(res, RadiationResult):
            force += np.array(list(res.force.values()))[None, :] * rao.sel(radiating_dof=res.radiating_dof).values[0, :, None]
        elif isinstance(res, DiffractionResult):
            idx = list(rao.coords["wave_direction"].values).index(res.wave_direction)
            force[0, idx, :] += np.array(list(res.force.values())) + np.array(list(froude_krylov_force(res).values()))

    return force[0, :, :] # forces[i_dir, influenced_dof]

def transformation_matrix(rao_k, rao_l):
    xk = rao_k.sel(radiating_dof="Roll")
    yk = rao_k.sel(radiating_dof="Pitch")
    zk = rao_k.sel(radiating_dof="Yaw")
    xl = rao_l.sel(radiating_dof="Roll")
    yl = rao_l.sel(radiating_dof="Pitch")
    zl = rao_l.sel(radiating_dof="Yaw")
    zero = np.zeros(xk.shape)

    H = (1/2) * (np.array([
                [-(yk * np.conjugate(yl) + zk * np.conjugate(zl)), zero, zero],
                [xk * np.conjugate(yl) + np.conjugate(xl) * yk, -(xk * np.conjugate(xl) + zk * np.conjugate(zl)), zero],
                [xk * np.conjugate(zl) + np.conjugate(xl) * zk, yk * np.conjugate(zl) + np.conjugate(yl) * zk, -(xk * np.conjugate(xl) + yk * np.conjugate(yl))]
            ],
            dtype=object)
        )

    H = np.reshape(H, (3, 3)) # H[xyz, xyz]
    return H

def skew_matrix(a):
    if a.shape[-1] != 3:
        raise ValueError("Last dimension should be of size 3.")

    results = np.zeros(a.shape[:-1] + (3,3), dtype=a.dtype)
    x, y, z = a[..., 0], a[..., 1], a[..., 2]
    results[..., 0, 1] = -z
    results[..., 0, 2] = y
    results[..., 1, 0] = z
    results[..., 1, 2] = -x
    results[..., 2, 0] = -y
    results[..., 2, 1] = x

    return results

def integrate_pressure_waterline(body, pressure):
        forces = {}
        normal = compute_faces_normals(body.mesh.vertices, body.mesh.faces_waterline)
        normal_waterline = normal/np.sqrt(1-normal[:, -1]**2)[:, None] # normal_waterline[i_edge_waterline, xyz]
        vertex_waterline = (
        body.mesh.vertices[body.mesh.edges_waterline[:, 0], :]
        + body.mesh.vertices[body.mesh.edges_waterline[:, 1], :]
    ) / 2
        for dof_name in body.dofs:
            dof = body.dofs[dof_name].evaluate_motion_at_points(vertex_waterline)

            # Scalar product on each edge:
            normal_dof_amplitude_on_waterline = np.sum(dof * normal_waterline, axis=1)
            forces[dof_name] = -np.sum(pressure * normal_dof_amplitude_on_waterline * body.mesh.length_edges_waterline)
        return forces
