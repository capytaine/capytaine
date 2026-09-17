"""Test related to the computation of mean drift force."""

import numpy as np
import xarray as xr
import pytest

import capytaine as cpt
from capytaine.io.xarray import problems_from_dataset, kochin_data_array
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force, near_field_mean_drift_force

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
    rao = cpt.post_pro.rao(dataset)
    print(rao)
    mdf = far_field_mean_drift_force(rao, dataset)
    force_analytical = dataset['g'].values * dataset['rho'].values * r * np.array([0.26, 0.7])
    assert np.allclose(mdf.sel(wave_direction_k=0, wave_direction_l=0)['drift_force_surge'], force_analytical, rtol=2e-1)
    assert "wave_direction_k" in mdf.dims
    assert "wave_direction_l" in mdf.dims
    assert tuple(mdf.sizes.values()) == (2, 2, 2)
    assert all(np.allclose(mdf.isel(wave_direction_k=0, wave_direction_l=1)[var].data,
                       mdf.isel(wave_direction_k=1, wave_direction_l=0).conj()[var].data) for var in mdf.data_vars)


def test_near_field_mean_drift_force():
    r = 1
    mesh = cpt.mesh_sphere(radius=r).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    wave_direction = [0, np.pi/4]
    k = np.array([0.92, 1.05])
    test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
        })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, results, solver)
    force_analytical = dataset['g'].values * dataset['rho'].values * r * np.array([0.26, 0.7])
    assert np.allclose(mdf.sel(wave_direction_k=0, wave_direction_l=0, influenced_dof='Surge'), force_analytical, rtol=2e-1)
    assert "wave_direction_k" in mdf.dims
    assert "wave_direction_l" in mdf.dims
    assert mdf.shape == (2, 2, 2, 6)  # (wavenumber, dir, dir, dofs)
    assert np.allclose(mdf.isel(wave_direction_k=0, wave_direction_l=1).values,
                        np.conj(mdf.isel(wave_direction_k=1, wave_direction_l=0).values))


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
        rao = cpt.post_pro.rao(dataset)
        force.append(far_field_mean_drift_force(rao, dataset)/r)

    assert np.allclose(force[0]['drift_force_surge'], force[1]['drift_force_surge'])


def test_scale_near_field_mean_drift_force():
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
        pbs = problems_from_dataset(test_matrix, body)
        results = solver.solve_all(pbs)
        dataset = cpt.assemble_dataset(results)
        rao = cpt.post_pro.rao(dataset)
        force.append(near_field_mean_drift_force(rao, results, solver)/r)

    assert np.allclose(force[0][...,0], force[1][...,0])


def test_cylinder_mean_drift_force():
    mesh = cpt.mesh_vertical_cylinder(length=1., resolution=(4,12,8)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    wave_direction = 27*np.pi/180
    omega = [3.]
    theta = np.linspace(-0.5, 2*np.pi, 20)
    test_matrix = xr.Dataset(coords={
                'omega': omega, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys()), 'theta': theta,
            })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    data_kochin = kochin_data_array(results, theta)
    dataset = cpt.assemble_dataset(results)
    dataset.update(data_kochin)
    rao = cpt.post_pro.rao(dataset)
    mdf_nf = near_field_mean_drift_force(rao, results, solver)/1e3
    mdf_ff = far_field_mean_drift_force(rao, dataset)/1e3
    target_fx = 3.22
    target_fy = 1.67
    target_fz = 12.49
    target_mx = -0.93
    target_my = 1.74
    target_mz = 0.

    assert np.isclose(mdf_nf[...,0], target_fx, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_ff['drift_force_surge'], target_fx, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_nf[...,1], target_fy, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_ff['drift_force_sway'], target_fy, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_nf[...,2], target_fz, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_nf[...,3], target_mx, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_nf[...,4], target_my, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_nf[...,5], target_mz, atol=1e-2, rtol=5e-1)
    assert np.isclose(mdf_ff['drift_force_yaw'], target_mz, atol=1e-2, rtol=5e-1)


def test_caisson():
    mesh = cpt.mesh_parallelepiped(size=(90,90,80)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    period = [14.20]
    theta = np.linspace(-0.5, 2*np.pi, 20)
    test_matrix = xr.Dataset(coords={
                'period': period, 'wave_direction': 0, 'radiating_dof': list(body.dofs.keys()), 'theta': theta,
            })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    data_kochin = kochin_data_array(results, theta)
    dataset = cpt.assemble_dataset(results)
    dataset.update(data_kochin)
    rao = cpt.post_pro.rao(dataset)
    mdf_nf = near_field_mean_drift_force(rao, results, solver)
    mdf_ff = far_field_mean_drift_force(rao, dataset)

    target_fx = 233204.32
    target_fz = 78824.94
    target_my = 10762390.49

    assert np.isclose(mdf_nf[..., 0], target_fx, rtol=4e-1)
    assert np.isclose(mdf_ff['drift_force_surge'], target_fx, rtol=4e-1)
    assert np.isclose(mdf_nf[..., 2], target_fz, rtol=4e-1)
    assert np.isclose(mdf_nf[..., 4], target_my, rtol=4e-1)


def test_symmetry_mean_drift_force():
    mesh = cpt.mesh_parallelepiped().immersed_part()
    mesh_sym = cpt.mesh_parallelepiped(reflection_symmetry=True).immersed_part()

    wave_direction = np.pi/4
    theta = np.linspace(-0.5, 2*np.pi, 20)
    k = np.array([2.5])
    solver = cpt.BEMSolver()

    mdf_ff = []
    mdf_nf = []
    for m in [mesh, mesh_sym]:
        body = cpt.FloatingBody(mesh=m, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
        body.inertia_matrix = body.compute_rigid_body_inertia()
        body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
        test_matrix = xr.Dataset(coords={
                'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
            })
        pbs = problems_from_dataset(test_matrix, body)
        results = solver.solve_all(pbs)
        data_kochin = kochin_data_array(results, theta)
        dataset = cpt.assemble_dataset(results)
        dataset.update(data_kochin)
        rao = cpt.post_pro.rao(dataset)
        mdf_nf.append(near_field_mean_drift_force(rao, results, solver))
        mdf_ff.append(far_field_mean_drift_force(rao, dataset))

    assert np.allclose(mdf_nf[0][...,0], mdf_nf[1][...,0])
    assert all(np.allclose(mdf_ff[0][var].data, mdf_ff[1][var].data) for var in mdf_ff[0].data_vars)


def test_period_equivalent_omega_mean_drift_force():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    solver = cpt.BEMSolver()
    wave_direction = [np.pi/3]
    period = np.array([1.6, 1.9, 2.4])

    test_matrix = xr.Dataset(coords={
            'period': period, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
        })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf_period = near_field_mean_drift_force(rao, results, solver)

    omega = 2*np.pi/period
    test_matrix = xr.Dataset(coords={
                'omega': omega, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
            })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf_omega = near_field_mean_drift_force(rao, results, solver)

    assert np.allclose(np.flip(mdf_period.values, axis=0), mdf_omega)
