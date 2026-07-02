"""Test related to the computation of mean drift force."""

import numpy as np
import xarray as xr
import pytest

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force

def test_far_field_mean_drift_force():
    r = 1
    mesh = cpt.mesh_sphere(radius=r).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    wave_direction = [0, np.pi/4]
    theta = np.linspace(-0.5, 2*np.pi, 20)
    k = np.array([0.92, 1.05])
    test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
        })
    dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True)
    X = cpt.post_pro.rao(dataset)
    res = far_field_mean_drift_force(X, dataset)['drift_force_surge'].squeeze()
    force_analytical = dataset['g'].values * dataset['rho'].values * r * np.array([0.26, 0.7])
    assert np.allclose(res.sel(wave_direction_k=0, wave_direction_l=0), force_analytical, rtol=2e-1)
    assert "wave_direction_k" in res.dims
    assert "wave_direction_l" in res.dims
    assert res.shape == (2, 2, 2)
    assert np.allclose(res.isel(wave_direction_k=0, wave_direction_l=1).values,
                       np.conj(res.isel(wave_direction_k=1, wave_direction_l=0).values))


def test_scale_far_field_mean_drift_force():
    radius = [1,5]
    force = []
    for r in radius:
        mesh = cpt.mesh_sphere(radius=r).immersed_part()
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
        body.inertia_matrix = body.compute_rigid_body_inertia()
        body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
        solver = cpt.BEMSolver()
        wave_direction = 0
        theta = np.linspace(-0.5, 2*np.pi, 20)
        k = np.linspace(0.1,1,5)/r
        test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
        })
        dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True)
        X = cpt.post_pro.rao(dataset)
        force.append(far_field_mean_drift_force(X, dataset)/r)

    assert np.allclose(force[0]['drift_force_surge'].values, force[1]['drift_force_surge'].values)
