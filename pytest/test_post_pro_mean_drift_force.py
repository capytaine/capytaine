"""Test related to the computation of mean drift force."""

import numpy as np
import xarray as xr
import pytest

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force, _merge_far_field_component

def test_far_field_mean_drift_force():
    r = 1
    mesh = cpt.mesh_sphere(radius=r).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    wave_direction = 0
    theta = np.linspace(0, 2*np.pi, 20)
    k = np.array([0.92, 1.05])
    test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
        })
    dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True)
    X = cpt.post_pro.rao(dataset)
    factor = dataset['g']*dataset['rho']*r
    res = far_field_mean_drift_force(X, dataset)['drift_force_surge'].squeeze()
    force_analytical = np.array([0.26, 0.7])
    assert np.allclose(res/factor, force_analytical, atol=1e-1)


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
        theta = np.linspace(0, 2*np.pi, 20)
        k = np.linspace(0.1,1,5)/r
        test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
        })
        dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True)
        X = cpt.post_pro.rao(dataset)
        force.append(far_field_mean_drift_force(X, dataset)/r)

    assert np.allclose(force[0]['drift_force_surge'].values, force[1]['drift_force_surge'].values)


def test_merge_far_field_component():
    """Test that _merge_far_field_component correctly merges the three force components."""
    # Create a simple mock dataset with the expected structure
    omega = np.array([1.0, 2.0])
    wave_direction = np.array([0.0, np.pi/4])

    surge_data = np.array([[1.0, 2.0], [3.0, 4.0]])
    sway_data = np.array([[5.0, 6.0], [7.0, 8.0]])
    yaw_data = np.array([[9.0, 10.0], [11.0, 12.0]])

    dataset = xr.Dataset({
        'drift_force_surge': xr.DataArray(surge_data, coords={'omega': omega, 'wave_direction': wave_direction}, dims=['omega', 'wave_direction']),
        'drift_force_sway': xr.DataArray(sway_data, coords={'omega': omega, 'wave_direction': wave_direction}, dims=['omega', 'wave_direction']),
        'drift_force_yaw': xr.DataArray(yaw_data, coords={'omega': omega, 'wave_direction': wave_direction}, dims=['omega', 'wave_direction']),
    })

    # Merge the components
    merged = _merge_far_field_component(dataset)

    # Check that the merged dataset has only one variable
    assert list(merged.data_vars) == ['drift_force']

    # Check that the influenced_dof dimension exists with correct coordinates
    assert 'influenced_dof' in merged['drift_force'].dims
    assert list(merged['drift_force'].influenced_dof.values) == ['Surge', 'Sway', 'Yaw']

    # Check that the values are correct
    assert np.allclose(merged['drift_force'].sel(influenced_dof='Surge').values, surge_data)
    assert np.allclose(merged['drift_force'].sel(influenced_dof='Sway').values, sway_data)
    assert np.allclose(merged['drift_force'].sel(influenced_dof='Yaw').values, yaw_data)
