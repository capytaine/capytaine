"""Computation of the mean drift force."""

# Copyright (C) 2026 the Capytaine developpers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import numpy as np
import xarray as xr
from collections import Counter

from capytaine.bem.airy_waves import (
    airy_waves_free_surface_elevation,
    airy_waves_velocity,
    froude_krylov_force,
)
from capytaine.bodies.dofs import AbstractDof
from capytaine.meshes.geometry import compute_faces_normals
from capytaine.bem.problems_and_results import DiffractionResult, RadiationResult


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
    k = dataset['wavenumber']
    h = dataset['water_depth']
    beta = dataset['wave_direction']
    rho = dataset['rho']
    g = dataset['g']
    H_diff = dataset['kochin_diffraction']
    H_rad = dataset['kochin_radiation'] 

    H_rad_tot = sum(H_rad.sel(radiating_dof=d)*X.sel(radiating_dof=d) for d in X.radiating_dof) 
    H_tot = np.exp(1j*np.pi/2)*(H_diff + H_rad_tot)
    H_beta = H_tot.interp(theta=beta)  

    H_derivative = H_tot.differentiate("theta") 
    H_derivative_beta = H_derivative.interp(theta=beta)

    coef1 = 2*np.pi*rho*omega
    if h == np.inf:
        coef2 = 2*np.pi*rho*k**2
    else:
        k0 = omega**2/g 
        coef2 = 2*np.pi*rho*k*(k0*h)**2 / (h*((k*h)**2 - (k0*h)**2 + k0*h))

    dims = [d for d in X.dims if d != 'radiating_dof']
    coords = {c: X.coords[c] for c in X.coords if c != 'radiating_dof'}

    base = np.abs(H_tot)**2
    
    Fx = xr.DataArray(
        data=-coef1*np.cos(beta)*np.imag(H_beta) - coef2*(base * np.cos(H_tot.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_surge'
        )
    
    Fy = xr.DataArray(
        data=-coef1*np.sin(beta)*np.imag(H_beta) - coef2*(base * np.sin(H_tot.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_sway'
        )
    
    Mz = xr.DataArray(
        data=coef1/k*np.real(H_derivative_beta) - coef2/k*np.imag((np.conjugate(H_tot) * H_derivative).integrate("theta")),
        coords=coords,
        dims=dims,
        name='drift_force_yaw'
        )
    
    return xr.Dataset({Fx.name: Fx, Fy.name: Fy, Mz.name: Mz})


def near_field_mean_drift_force(rao, results, solver):
    body = results[0].body
    mesh = body.mesh
    rho = results[0].rho
    g = results[0].g

    nb_dir = rao.sizes["wave_direction"]
    zero_block = np.zeros((nb_dir, 3, 3))

    z1 = motion_order1(mesh.faces_centers, rao)[:, :, -1] # z1[i_dir, i_face]
    forces_order1 = hydrodynamics_forces_order1(results, rao) + hydrostatics_forces(body, rao, rho, g, z1)
    rotation = rao.sel(radiating_dof=["Roll", "Pitch", "Yaw"]).values[0, :, :] # rotation[i_dir, xyz]
    rotation_matrix = np.block([[skew_matrix(rotation), zero_block], [zero_block, skew_matrix(rotation)]])
    rotation_forces_order1 = (1/2) * np.real(np.matvec(rotation_matrix, np.conjugate(forces_order1)))

    z0 = np.tile(mesh.faces_centers[:, -1], (nb_dir, 1))
    forces_order0 = hydrostatics_forces(body, rao, rho, g, z0)
    H = np.block([[transformation_matrix(rao), zero_block], [transformation_matrix(rao), zero_block]])
    hydrostatics_order2 = (1/2) * np.matvec(H, forces_order0)

    translation = rao.sel(radiating_dof=["Surge", "Sway", "Heave"]).values[0, :, :] # translation[i_dir, xyz]
    translation_matrix = np.block([[zero_block, zero_block], [skew_matrix(translation), zero_block]])
    rotation_forces_order0 = np.matvec(rotation_matrix, forces_order0)
    translation_moment = (1/2) * np.real(np.matvec(translation_matrix, np.conjugate(forces_order1 + rotation_forces_order0)))

    all_terms = rotation_forces_order1 + hydrostatics_order2 + translation_moment
    pressure_field = pressure_field_order2(mesh, solver, results, rao, g)
    field_waterline = pertubated_part(mesh, solver, results, rao, g)

    for idx, beta in enumerate(rao.coords["wave_direction"].values):
        pressure_hull = body.integrate_pressure(pressure_field.sel(wave_direction=beta))
        pressure_waterline = integrate_pressure_waterline(body, field_waterline.sel(wave_direction=beta))
        all_terms[idx, :] += rho * (np.array(list(pressure_hull.values())) + np.array(list(pressure_waterline.values())))

    coords = {
        "wave_direction": rao.coords["wave_direction"],
        "influenced_dof": list(pressure_hull.keys()),
    }

    F = xr.DataArray(data=all_terms, coords=coords)
    return F

def motion_order1(points, rao):
    translation = rao.sel(radiating_dof=["Surge", "Sway", "Heave"]).values
    rotation = rao.sel(radiating_dof=["Roll", "Pitch", "Yaw"]).values
    results = translation[0, :, None, :] + np.cross(rotation[0, :, None, :], points[None, :, :])
    # results[i_dir, i_face, xyz] = translation[i_dir, xyz] + cross(rotation[i_dir, xyz], point[i_face, xyz])
    return results

def pressure_field_order2(mesh, solver, results, rao, g):
    gradient_potential = total_potential_gradient(solver, mesh, results, rao)
    motion = motion_order1(mesh.faces_centers, rao)
    product = np.real(np.sum(motion * np.conjugate(-1j * rao.coords["omega"].values * gradient_potential), axis=2))
    gradient_potential_square = np.sum(np.abs(gradient_potential) ** 2, axis=2)
    hydrostatic = np.matvec(transformation_matrix(rao)[:, None, :, :], mesh.faces_centers[None, :, :])
    # hydrostatic[i_dir, i_face, xyz] = matvec(transofrmation_matrix[i_dir, xyz, xyz], faces_centers[i_face, xyz])
    z_order2 = hydrostatic[:, :, -1] # z_order2[i_dir, i_face]
    return - (1/2) * (product + gradient_potential_square/2 + g * z_order2)

def pertubated_part(mesh, solver, results, rao, g):
    edges_waterline = mesh.edges_waterline
    vertices_middle_waterline = (
        mesh.vertices[edges_waterline[:, 0], :]
        + mesh.vertices[edges_waterline[:, 1], :]
    ) / 2
    free_surface_elevation = total_free_surface_elevation(solver, vertices_middle_waterline, results, rao)
    vertical_position = motion_order1(vertices_middle_waterline, rao)[:, :, -1] # vertical_postion[i_dir, i_face]
    return (1/4) * g * np.abs(free_surface_elevation - vertical_position) ** 2

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
            force += list(res.force.values()) * rao
        elif isinstance(res, DiffractionResult):
            force += np.array(list(res.force.values())) + np.array(list(froude_krylov_force(res).values()))

    return force[0, :, :].values # forces[i_dir, influenced_dof]

def transformation_matrix(rao):
    x = rao.sel(radiating_dof="Roll")
    y = rao.sel(radiating_dof="Pitch")
    z = rao.sel(radiating_dof="Yaw")
    zero = np.zeros(x.shape)

    H = (1/2) * (np.array([
                [-(np.abs(y) ** 2 + np.abs(z) ** 2), zero, zero],  
                [2 * np.real(x * np.conjugate(y)), -(np.abs(x) ** 2 + np.abs(z) ** 2), zero],
                [2 * np.real(x * np.conjugate(z)), 2 * np.real(y * np.conjugate(z)), -(np.abs(x) ** 2 + np.abs(y) ** 2)]
            ],
            dtype=object)
        )
    
    H = np.reshape(H, (rao.sizes['wave_direction'], 3, 3)) # H[i_dir, xyz, xyz]
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
        normal[:, -1] = 0
        normal_waterline = normal/np.linalg.norm(normal[:, :2], ord=2, axis=1, keepdims=True)
        vertex_waterline = np.sort(body.mesh.edges_waterline[:,0])
        for dof_name in body.dofs:
            if isinstance(body.dofs[dof_name], AbstractDof):
                dof = body.dofs[dof_name].evaluate_motion(body.mesh)
            else:
                dof = body.dofs[dof_name]
            # Scalar product on each edge:
            normal_dof_amplitude_on_waterline = np.sum(dof[vertex_waterline] * normal_waterline, axis=1)
            forces[dof_name] = -np.sum(pressure * normal_dof_amplitude_on_waterline * body.mesh.length_edges_waterline)
        return forces