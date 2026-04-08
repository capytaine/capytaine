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
    mesh = results[0].body.mesh
    rho = results[0].rho
    g = results[0].g

    total_free_surf_elev = free_surface_elevation_field(mesh, solver, results, rao)
    normal = compute_faces_normals(mesh.vertices, mesh.faces_water_line)
    normal[:, -1] = 0
    normal_waterline = normal/np.linalg.norm(normal[:, :2], ord=2, axis=1, keepdims=True)
    data = total_free_surf_elev.values[:, :, None] * normal_waterline[None, :, :]
    # data[i_dir, i_face_waterline, xyz] = total_free_surf_elev[i_dir, i_face_waterline] * normal_waterline[i_face_waterline, xyz]
    waterline_integral = -(rho * g / 4) * np.sum(data * mesh.length_edges_water_line[None, :, None], axis=1) # sum over the waterline edges 

    forces_order1 = hydrodynamics_forces_order1(results, rao) + hydrostatics_forces_order1(mesh, rao, rho, g)
    rotation = rao.sel(radiating_dof=["Roll", "Pitch", "Yaw"]).values[0, :, :] # rotation[i_dir, xyz]
    rotation_forces_order1 = (1/2) * np.real(np.cross(rotation, np.conjugate(forces_order1), axis=1))

    z0 = mesh.faces_centers[:, -1]
    H = transformation_matrix(rao)
    n = mesh.faces_normals
    data = z0[None, :, None] * np.matvec(H[:, None, :, :], n[None, :, :])
    # data[i_dir, i_face, xyz] = z0[i_face] * matvec(H[i_dir, xyz, xyz], n[i_face, xyz])
    hydrostatics_order2 = -(1/4) * rho * g * np.sum(data * mesh.faces_areas[None, : , None], axis=1)

    all_terms = waterline_integral + rotation_forces_order1 + hydrostatics_order2
    pressure_field = pressure_field_order2(mesh, solver, results, rao, g)

    for idx, beta in enumerate(rao.coords["wave_direction"].values):
        pressure_hull = results[0].body.integrate_pressure(pressure_field.sel(wave_direction=beta))
        all_terms[idx, :] += rho * np.array(list(pressure_hull.values()))[:3]

    coords = {
        "wave_direction": rao.coords["wave_direction"],
        "influenced_dof": list(pressure_hull.keys())[:3],
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
    product = np.real(np.sum(motion* np.conjugate(-1j * rao.coords["omega"].values * gradient_potential),axis=2))
    gradient_potential_square = np.sum(np.abs(gradient_potential) ** 2, axis=2)
    hydrostatic = np.matvec(transformation_matrix(rao)[:, None, :, :], mesh.faces_centers[None, :, :])
    # hydrostatic[i_dir, i_face, xyz] = matvec(transofrmation_matrix[i_dir, xyz, xyz], faces_centers[i_face, xyz])
    z_order2 = hydrostatic[:, :, -1] # z_order2[i_dir, i_face]
    return (
        -gradient_potential_square / 4
        - product / 2
        - g * z_order2 / 2
    )

def free_surface_elevation_field(mesh, solver, results, X):
    edges_waterline = mesh.edges_waterline
    vertices_middle_waterline = (
        mesh.vertices[edges_waterline[:, 0], :]
        + mesh.vertices[edges_waterline[:, 1], :]
    ) / 2
    free_surface_elevation = total_free_surface_elevation(solver, vertices_middle_waterline, results, X)
    vertical_position = motion_order1(vertices_middle_waterline, X)[:, :, -1] # vertical_postion[i_dir, i_face]
    return np.abs(free_surface_elevation - vertical_position) ** 2

def total_potential_gradient(solver, mesh, results, X):
    potential_gradient = xr.DataArray(
        data=np.zeros((X.sizes["wave_direction"], mesh.nb_faces, 3), dtype=complex),
        coords={
            "wave_direction": X.coords["wave_direction"],
            "mesh_face": np.arange(mesh.nb_faces),
            "xyz": np.arange(3),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            potential_gradient.loc[{"wave_direction": res.wave_direction}] += airy_waves_velocity(mesh, res) # incident
            potential_gradient.loc[{"wave_direction": res.wave_direction}] += solver._compute_potential_gradient(mesh, res) # diffracted

        elif isinstance(res, RadiationResult):
            X_rad = X.sel({"radiating_dof": res.radiating_dof, "wave_direction": res.wave_direction,}).values
            potential_gradient.loc[{"wave_direction": res.wave_direction}] += X_rad * solver._compute_potential_gradient(mesh, res) # radiated

    return potential_gradient

def total_free_surface_elevation(solver, vertices_middle, results, X):
    nb_vertices_waterline = np.shape(vertices_middle)[0]
    free_surface_elevation = xr.DataArray(
        data=np.zeros((X.sizes["wave_direction"], nb_vertices_waterline), dtype=complex),
        coords={
            "wave_direction": X.coords["wave_direction"],
            "vertices_waterline": np.arange(nb_vertices_waterline),
        },
    )
    for res in results:
        if isinstance(res, DiffractionResult):
            free_surface_elevation.loc[{"wave_direction": res.wave_direction}] += airy_waves_free_surface_elevation(vertices_middle, res) # incident
            free_surface_elevation.loc[{"wave_direction": res.wave_direction}] += solver.compute_free_surface_elevation(vertices_middle, res) # diffracted

        elif isinstance(res, RadiationResult):
            X_rad = X.sel({"radiating_dof": res.radiating_dof, "wave_direction": res.wave_direction}).values
            free_surface_elevation.loc[{"wave_direction": res.wave_direction}] += X_rad * solver.compute_free_surface_elevation(vertices_middle, res) # radiated

    return free_surface_elevation

def hydrostatics_forces_order1(mesh, rao, rho, g):
    z1 = motion_order1(mesh.faces_centers, rao)[:, :, -1] 
    data = z1[:, :, None] * mesh.faces_normals[None, :, :]
    # data[i_dir, i_face, xyz] = z1[i_dir, i_face] * faces_normals[i_face, xyz]
    results = rho * g * np.sum(data * mesh.faces_areas[None, :, None], axis=1)
    return results

def hydrodynamics_forces_order1(results, rao):
    force = np.zeros(list(rao.sizes.values()), dtype=complex)

    for res in results:
        if isinstance(res, RadiationResult):
            force += list(res.force.values()) * rao
        elif isinstance(res, DiffractionResult):
            force += np.array(list(res.force.values())) + np.array(list(froude_krylov_force(res).values()))

    return force.sel(radiating_dof=["Surge", "Sway", "Heave"])[0, :, :] # forces[i_dir, radiating_dof in ["Surge", "Sway", "Heave"]]

def transformation_matrix(rao):
    x = rao.sel(radiating_dof="Roll")
    y = rao.sel(radiating_dof="Pitch")
    z = rao.sel(radiating_dof="Yaw")
    zero = np.zeros(x.shape)

    H = 0.5*(np.array([
                [-(np.abs(y) ** 2 + np.abs(z) ** 2), zero, zero],  
                [2 * np.real(x * np.conjugate(y)), -(np.abs(x) ** 2 + np.abs(z) ** 2), zero],
                [2 * np.real(x * np.conjugate(z)), 2 * np.real(y * np.conjugate(z)), -(np.abs(x) ** 2 + np.abs(y) ** 2)]
            ],
            dtype=object)
        )
    
    H = np.reshape(H, (rao.sizes['wave_direction'], 3, 3)) # H[i_dir, xyz, xyz]
    return H 
