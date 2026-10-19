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
"""Test related to the computation of mean drift force."""

import numpy as np
import xarray as xr
import pytest

import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.io.xarray import problems_from_dataset, kochin_data_array
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force, near_field_mean_drift_force, _near_field_hydrodynamics


def test_far_field_mean_drift_force():
    r = 1
    mesh = cpt.mesh_sphere(radius=r).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    solver = helpers.solver()
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
    solver = helpers.solver()
    wave_direction = [0, np.pi/4]
    k = np.array([0.92, 1.05])
    test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
        })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, results, solver)['near_field_mean_drift_force']
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
        solver = helpers.solver()
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
        solver = helpers.solver()
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
        force.append(near_field_mean_drift_force(rao, results, solver)['near_field_mean_drift_force']/r)

    assert np.allclose(force[0][...,0], force[1][...,0])


def test_cylinder_mean_drift_force():
    mesh = cpt.mesh_vertical_cylinder(length=1., resolution=(4,12,8)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    solver = helpers.solver()
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
    mdf_nf = near_field_mean_drift_force(rao, results, solver)['near_field_mean_drift_force'].values/1e3
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
    solver = helpers.solver()
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
    mdf_nf = near_field_mean_drift_force(rao, results, solver)['near_field_mean_drift_force'].values
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
    solver = helpers.solver()

    mdf_ff = []
    mdf_nf = []
    for m in [mesh, mesh_sym]:
        body = cpt.FloatingBody(mesh=m, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
        test_matrix = xr.Dataset(coords={
                'wavenumber': k, 'wave_direction': wave_direction, 'theta': theta, 'radiating_dof': list(body.dofs.keys())
            })
        pbs = problems_from_dataset(test_matrix, body)
        results = solver.solve_all(pbs)
        data_kochin = kochin_data_array(results, theta)
        dataset = cpt.assemble_dataset(results)
        dataset.update(data_kochin)
        rao = cpt.post_pro.rao(dataset)
        mdf_nf.append(near_field_mean_drift_force(rao, results, solver)['near_field_mean_drift_force'])
        mdf_ff.append(far_field_mean_drift_force(rao, dataset))

    assert np.allclose(mdf_nf[0][...,0], mdf_nf[1][...,0])
    assert all(np.allclose(mdf_ff[0][var].data, mdf_ff[1][var].data) for var in mdf_ff[0].data_vars)


def test_period_equivalent_omega_mean_drift_force():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    solver = helpers.solver()
    wave_direction = [np.pi/3]
    period = np.array([1.6, 1.9, 2.4])

    test_matrix = xr.Dataset(coords={
            'period': period, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
        })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf_period = near_field_mean_drift_force(rao, results, solver)["near_field_mean_drift_force"]

    omega = 2*np.pi/period
    test_matrix = xr.Dataset(coords={
                'omega': omega, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
            })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf_omega = near_field_mean_drift_force(rao, results, solver)["near_field_mean_drift_force"]

    assert np.allclose(np.flip(mdf_period.values, axis=0), mdf_omega)


def test_near_field_hydrodynamics_of_generic_dofs():
    # Smoke test for dofs that are not rigid body motions (no reference value).
    bending = cpt.CustomDof(motion=lambda x: np.array([0.0, 0.0, x[0]**2]),
                            gradient_of_motion=lambda x: np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2*x[0], 0.0, 0.0]]))
    body = cpt.FloatingBody(mesh=helpers.small_sphere_mesh(),
                            dofs={**cpt.rigid_body_dofs(only=["Heave", "Pitch"], rotation_center=(0, 0, 0)), "Bending": bending})
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0, np.pi/2], 'radiating_dof': list(body.dofs)})
    results = solver.solve_all(problems_from_dataset(test_matrix, body))
    rng = np.random.default_rng(3)
    rao = xr.DataArray(rng.normal(size=(1, 2, 3)) + 1j * rng.normal(size=(1, 2, 3)),
                       dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0], "wave_direction": [0.0, np.pi/2], "radiating_dof": list(body.dofs)})
    force, pressure_hull, pressure_waterline, relative_elevation = _near_field_hydrodynamics(rao, results, solver)
    assert force.shape == (1, 2, 2, 3)
    assert np.all(np.isfinite(force))
    np.testing.assert_allclose(force, np.conj(np.swapaxes(force, 1, 2)), atol=1e-10)


################################################################################
#                          Mean drift force pressure                           #
################################################################################

def test_mean_drift_force_pressure():
    body = helpers.small_sphere_body("rigid", lid=True)
    solver = helpers.solver()
    wave_direction = [0]
    k = np.array([0.92])
    test_matrix = xr.Dataset(coords={
            'wavenumber': k, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
        })
    pbs = problems_from_dataset(test_matrix, body)
    results = solver.solve_all(pbs)
    dataset = cpt.assemble_dataset(results)
    rao = cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, results, solver, output_pressure=True)
    assert mdf["pressure_hull_mean_drift"].shape == (1, 1, 1, body.mesh.nb_faces)


def test_mean_drift_force_from_exported_pressures_of_fixed_body():
    # Without body motion, the near-field force is the integral of the second-order
    # pressure on the hull, plus the waterline term, for each pair of wave directions.
    body = helpers.small_sphere_body("rigid", lid=True)
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
            'wavenumber': [0.92], 'wave_direction': [0.0, np.pi/3], 'radiating_dof': list(body.dofs.keys())
        })
    results = solver.solve_all(problems_from_dataset(test_matrix, body))
    rao = 0.0 * cpt.post_pro.rao(cpt.assemble_dataset(results))
    mdf = near_field_mean_drift_force(rao, results, solver, output_pressure=True).isel(wavenumber=0)

    p_wl = mdf["pressure_waterline_total"].transpose("wave_direction", "waterline_segment").values
    h_wl = mdf["relative_elevation_waterline"].transpose("wave_direction", "waterline_segment").values
    for k in range(2):
        for l in range(2):
            hull = body.integrate_pressure(mdf["pressure_hull_mean_drift"].isel(wave_direction_k=k, wave_direction_l=l).values)
            waterline = body.integrate_pressure_on_waterline(
                    (p_wl[k] * np.conjugate(h_wl[l]) + np.conjugate(p_wl[l]) * h_wl[k])/4
                    )
            force = mdf["near_field_mean_drift_force"].isel(wave_direction_k=k, wave_direction_l=l).values
            np.testing.assert_allclose(
                    np.array(list(hull.values())) + np.array(list(waterline.values())),
                    force, rtol=1e-10, atol=1e-10)
