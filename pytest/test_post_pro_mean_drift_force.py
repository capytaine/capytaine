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
from capytaine.io.xarray import separate_complex_values, merge_complex_values, save_dataset_as_netcdf
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force, near_field_mean_drift_force, _sources_on_mesh_including_lid
from capytaine.tools.perturbation_expansion import DifferenceFrequencyExpansion as Expansion


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
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, dataset, body, solver=solver)['near_field_mean_drift_force']
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
        dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
        rao = cpt.post_pro.rao(dataset)
        force.append(near_field_mean_drift_force(rao, dataset, body, solver=solver)['near_field_mean_drift_force']/r)

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
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf_nf = near_field_mean_drift_force(rao, dataset, body, solver=solver)['near_field_mean_drift_force'].values/1e3
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
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf_nf = near_field_mean_drift_force(rao, dataset, body, solver=solver)['near_field_mean_drift_force'].values
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
        dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
        rao = cpt.post_pro.rao(dataset)
        mdf_nf.append(near_field_mean_drift_force(rao, dataset, body, solver=solver)['near_field_mean_drift_force'])
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
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf_period = near_field_mean_drift_force(rao, dataset, body, solver=solver)["near_field_mean_drift_force"]

    omega = 2*np.pi/period
    test_matrix = xr.Dataset(coords={
                'omega': omega, 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs.keys())
            })
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf_omega = near_field_mean_drift_force(rao, dataset, body, solver=solver)["near_field_mean_drift_force"]

    assert np.allclose(np.flip(mdf_period.values, axis=0), mdf_omega)


def test_near_field_mean_drift_force_of_generic_dofs():
    # Smoke test for dofs that are not rigid body motions (no reference value).
    bending = cpt.CustomDof(motion=lambda x: np.array([0.0, 0.0, x[0]**2]),
                            gradient_of_motion=lambda x: np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2*x[0], 0.0, 0.0]]))
    body = cpt.FloatingBody(mesh=helpers.small_sphere_mesh(),
                            dofs={**cpt.rigid_body_dofs(only=["Heave", "Pitch"], rotation_center=(0, 0, 0)), "Bending": bending})
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0, np.pi/2], 'radiating_dof': list(body.dofs)})
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    rng = np.random.default_rng(3)
    rao = xr.DataArray(rng.normal(size=(1, 2, 3)) + 1j * rng.normal(size=(1, 2, 3)),
                       dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0], "wave_direction": [0.0, np.pi/2], "radiating_dof": list(body.dofs)})
    force = near_field_mean_drift_force(rao, dataset, body, solver=solver)["near_field_mean_drift_force"].values
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
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, dataset, body, solver=solver, output_pressure=True)
    assert mdf["pressure_hull_mean_drift"].shape == (1, 1, 1, body.mesh.nb_faces)


def test_mean_drift_force_from_exported_pressures_of_fixed_body():
    # Without body motion, the near-field force is the integral of the second-order
    # pressure on the hull, plus the waterline term, for each pair of wave directions.
    body = helpers.small_sphere_body("rigid", lid=True)
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
            'wavenumber': [0.92], 'wave_direction': [0.0, np.pi/3], 'radiating_dof': list(body.dofs.keys())
        })
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
    rao = 0.0 * cpt.post_pro.rao(dataset)
    mdf = near_field_mean_drift_force(rao, dataset, body, solver=solver, output_pressure=True).isel(wavenumber=0)

    p_wl = mdf["pressure_waterline_total"].transpose("wave_direction", "waterline_segment").values
    h_wl = mdf["relative_elevation_waterline"].transpose("wave_direction", "waterline_segment").values
    strip_width = h_wl / np.sqrt(1 - body.mesh.merged().normals_edges_waterline[:, 2]**2)
    waterline = body.integrate_pressure_on_waterline(Expansion.first_order(p_wl) * Expansion.first_order(strip_width)).coefficients[2]
    for k in range(2):
        for l in range(2):
            hull = body.integrate_pressure(mdf["pressure_hull_mean_drift"].isel(wave_direction_k=k, wave_direction_l=l).values)
            force = mdf["near_field_mean_drift_force"].isel(wave_direction_k=k, wave_direction_l=l).values
            np.testing.assert_allclose(np.array(list(hull.values())) + waterline[k, l], force, rtol=1e-10, atol=1e-10)


def test_near_field_mean_drift_force_of_submerged_body():
    # No waterline: the waterline term and the waterline outputs are empty.
    mesh = cpt.mesh_sphere(radius=0.5, center=(0.0, 0.0, -1.5), resolution=(4, 4))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=["Heave", "Pitch"], rotation_center=(0.0, 0.0, -1.5)))
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0], 'radiating_dof': list(body.dofs)})
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    rao = xr.DataArray(np.full((1, 1, 2), 1.0 + 0.5j), dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0], "wave_direction": [0.0], "radiating_dof": list(body.dofs)})
    mdf = near_field_mean_drift_force(rao, dataset, body, solver=solver, output_pressure=True)
    assert np.all(np.isfinite(mdf["near_field_mean_drift_force"].values))
    assert mdf["pressure_waterline_total"].sizes["waterline_segment"] == 0


def test_near_field_mean_drift_force_of_array_dofs():
    mesh = helpers.small_sphere_mesh()
    body = cpt.FloatingBody(mesh=mesh, dofs={"Heave": np.tile([0.0, 0.0, 1.0], (mesh.nb_faces, 1))})
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0], 'radiating_dof': ["Heave"]})
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    rao = xr.DataArray(np.ones((1, 1, 1), dtype=complex), dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0], "wave_direction": [0.0], "radiating_dof": ["Heave"]})
    with pytest.raises(NotImplementedError):
        near_field_mean_drift_force(rao, dataset, body, solver=solver)


def dataset_and_rao_of_sphere_with_lid():
    body = helpers.small_sphere_body("rigid", lid=True)
    test_matrix = xr.Dataset(coords={'omega': [1.0, 1.5], 'wave_direction': [0.0, np.pi/3], 'radiating_dof': list(body.dofs)})
    dataset = helpers.solver().fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    rng = np.random.default_rng(7)
    rao = xr.DataArray(rng.normal(size=(2, 2, 6)) + 1j * rng.normal(size=(2, 2, 6)), dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0, 1.5], "wave_direction": [0.0, np.pi/3], "radiating_dof": list(body.dofs)})
    return body, dataset, rao


def test_near_field_mean_drift_force_after_netcdf_round_trip(tmp_path):
    body, dataset, rao = dataset_and_rao_of_sphere_with_lid()
    path = tmp_path / "dataset.nc"
    save_dataset_as_netcdf(path, separate_complex_values(dataset))
    loaded = merge_complex_values(xr.load_dataset(path))
    force = near_field_mean_drift_force(rao, dataset, body, solver=helpers.solver())["near_field_mean_drift_force"]
    loaded_force = near_field_mean_drift_force(rao, loaded, body, solver=helpers.solver())["near_field_mean_drift_force"]
    np.testing.assert_allclose(loaded_force.values, force.values, rtol=1e-12)


def test_near_field_mean_drift_force_without_sources():
    body, dataset, rao = dataset_and_rao_of_sphere_with_lid()
    with pytest.raises(ValueError, match="keep_details=True"):
        near_field_mean_drift_force(rao, dataset.drop_vars("sources_lid_radiation"), body, solver=helpers.solver())
    with pytest.raises(ValueError, match="direct method"):
        near_field_mean_drift_force(rao, dataset.assign_attrs(method="direct"), body, solver=helpers.solver())


def test_sources_of_dataset_on_interleaved_hull_and_lid():
    # For a symmetric mesh with a lid, the faces of the hull and of the lid alternate in the mesh including lid.
    body = helpers.small_sphere_body("rigid", plane_symmetry=True, lid=True)
    assert np.any(np.diff(body.hull_mask.astype(int)) > 0)  # Some lid faces are before some hull faces
    problems = [cpt.DiffractionProblem(body=body, omega=1.0, wave_direction=beta) for beta in [0.0, np.pi/3]]
    problems += [cpt.RadiationProblem(body=body, omega=1.0, radiating_dof=dof) for dof in body.dofs]
    results = helpers.solver().solve_all(problems, keep_details=True)
    dataset = cpt.assemble_dataset(results, hydrostatics=False).isel(omega=0)
    sources = _sources_on_mesh_including_lid(dataset, "diffraction", "wave_direction", body)
    np.testing.assert_allclose(sources, np.array([r.sources for r in results[:2]]), rtol=1e-12)
    sources = _sources_on_mesh_including_lid(dataset.sel(radiating_dof=list(body.dofs)), "radiation", "radiating_dof", body)
    np.testing.assert_allclose(sources, np.array([r.sources for r in results[2:]]), rtol=1e-12)


def test_near_field_mean_drift_force_of_symmetric_body_with_lid():
    # The faces of the hull and of the lid alternate in the mesh including lid (see above):
    # the self-influence term of the velocity on the hull should still be added on the right panels.
    body_sym = helpers.small_sphere_body("rigid", plane_symmetry=True, lid=True)
    body = cpt.FloatingBody(mesh=body_sym.mesh.merged(), lid_mesh=body_sym.lid_mesh.merged(),
                            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)), center_of_mass=(0, 0, 0))
    forces = []
    for b in [body_sym, body]:
        test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0, np.pi/3], 'radiating_dof': list(b.dofs)})
        dataset = helpers.solver().fill_dataset(test_matrix, b, keep_details=True, hydrostatics=False)
        rao = xr.DataArray(np.full((1, 2, 6), 0.1 + 0.2j), dims=["omega", "wave_direction", "radiating_dof"],
                           coords={"omega": [1.0], "wave_direction": [0.0, np.pi/3], "radiating_dof": list(b.dofs)})
        forces.append(near_field_mean_drift_force(rao, dataset, b, solver=helpers.solver())["near_field_mean_drift_force"])
    np.testing.assert_allclose(forces[0].values, forces[1].values, rtol=1e-6, atol=1e-6*np.abs(forces[1].values).max())


def test_near_field_mean_drift_force_with_inconsistent_inputs():
    body, dataset, rao = dataset_and_rao_of_sphere_with_lid()
    solver = helpers.solver()
    body_without_lid = cpt.FloatingBody(mesh=body.mesh, dofs=body.dofs, name=body.name)
    with pytest.raises(ValueError, match="lid"):
        near_field_mean_drift_force(rao, dataset, body_without_lid, solver=solver)
    with pytest.raises(ValueError, match="no common omega"):
        near_field_mean_drift_force(rao.assign_coords(omega=[1.1, 1.6]), dataset, body, solver=solver)
    with pytest.raises(ValueError, match="no common radiating_dof"):
        near_field_mean_drift_force(rao.isel(radiating_dof=slice(0, 3)), dataset.isel(radiating_dof=slice(3, None)), body, solver=solver)
    with pytest.raises(NotImplementedError, match="single value of rho"):
        near_field_mean_drift_force(rao, xr.concat([dataset, dataset.assign_coords(rho=1025.0)], dim="rho"), body, solver=solver)


def test_near_field_mean_drift_force_on_intersection_of_inputs(caplog):
    body, dataset, rao = dataset_and_rao_of_sphere_with_lid()
    solver = helpers.solver()
    with caplog.at_level("WARNING"):
        reference = near_field_mean_drift_force(rao, dataset, body, solver=solver)["near_field_mean_drift_force"]
    assert "common" not in caplog.text
    with caplog.at_level("WARNING"):
        force = near_field_mean_drift_force(rao.assign_coords(omega=[1.0, 1.6]), dataset, body, solver=solver)["near_field_mean_drift_force"]
    assert "same omega" in caplog.text
    np.testing.assert_allclose(force.values, reference.sel(omega=[1.0]).values, rtol=1e-12)
    caplog.clear()
    dofs = list(body.dofs)
    with caplog.at_level("WARNING"):
        force = near_field_mean_drift_force(rao.sel(radiating_dof=dofs[:4]), dataset.sel(radiating_dof=dofs[1:]), body, solver=solver)["near_field_mean_drift_force"]
    assert "same radiating_dof" in caplog.text
    assert list(force.coords["influenced_dof"].values) == dofs[1:4]
    expected = near_field_mean_drift_force(rao.sel(radiating_dof=dofs[1:4]), dataset, body.with_only_dofs(dofs[1:4]), solver=solver)["near_field_mean_drift_force"]
    np.testing.assert_allclose(force.values, expected.values, rtol=1e-12)


def test_near_field_mean_drift_force_with_rao_in_another_order():
    body, dataset, rao = dataset_and_rao_of_sphere_with_lid()
    force = near_field_mean_drift_force(rao, dataset, body, solver=helpers.solver())["near_field_mean_drift_force"]
    reordered_rao = rao.isel(omega=[1, 0], wave_direction=[1, 0], radiating_dof=[5, 3, 1, 0, 2, 4]).transpose("radiating_dof", "omega", "wave_direction")
    reordered_force = near_field_mean_drift_force(reordered_rao, dataset, body, solver=helpers.solver())["near_field_mean_drift_force"]
    xr.testing.assert_allclose(reordered_force.sel(omega=force.omega, wave_direction_k=force.wave_direction_k, wave_direction_l=force.wave_direction_l), force)


def near_field_mean_drift_force_with_imposed_motion(body, motion_amplitude, wave_direction):
    """Near-field QTF of `body` at omega = 1 for an imposed motion `motion_amplitude[i_dir, i_dof]`."""
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': wave_direction, 'radiating_dof': list(body.dofs)})
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    rao = xr.DataArray(motion_amplitude[None], dims=["omega", "wave_direction", "radiating_dof"],
                       coords={"omega": [1.0], "wave_direction": wave_direction, "radiating_dof": list(body.dofs)})
    return near_field_mean_drift_force(rao, dataset, body, solver=solver)["near_field_mean_drift_force"].isel(omega=0)


def test_near_field_mean_drift_force_of_multibody_of_single_body():
    body = cpt.FloatingBody(mesh=helpers.small_sphere_mesh(), dofs=cpt.rigid_body_dofs(rotation_center=(0.1, 0.0, -0.2)), name="a")
    rng = np.random.default_rng(5)
    motion_amplitude = rng.normal(size=(2, 6)) + 1j * rng.normal(size=(2, 6))
    wave_direction = [0.0, np.pi/3]
    force = near_field_mean_drift_force_with_imposed_motion(body, motion_amplitude, wave_direction)
    multibody_force = near_field_mean_drift_force_with_imposed_motion(cpt.Multibody([body]), motion_amplitude, wave_direction)
    assert list(multibody_force.influenced_dof.values) == [f"a__{dof}" for dof in body.dofs]
    np.testing.assert_allclose(multibody_force.values, force.values, rtol=1e-10, atol=1e-10)


def test_near_field_mean_drift_force_of_mirrored_bodies():
    # Two bodies mirror images of each other through the plane yOz, with waves along y:
    # the problem is symmetric, so the forces on the two bodies are mirror images of each other,
    # interactions and second-order motion of each body included.
    # Planar faces: the first-order interactions of meshes with non-planar faces (e.g. a sphere) are only symmetric up to ~1e-4.
    mesh = cpt.mesh_parallelepiped(size=(1.0, 1.2, 1.0), resolution=(3, 3, 3), center=(-2.0, 0.1, -0.5)).immersed_part()
    a = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(-2.0, 0.0, -0.2)), name="a")
    b = cpt.FloatingBody(mesh=mesh.mirrored("yOz"), dofs=cpt.rigid_body_dofs(rotation_center=(2.0, 0.0, -0.2)), name="b")
    # Under the mirror symmetry x -> -x, the translations transform as vectors and the rotations as pseudo-vectors.
    mirror = np.array([-1.0, 1.0, 1.0, 1.0, -1.0, -1.0])  # Surge, Sway, Heave, Roll, Pitch, Yaw
    rng = np.random.default_rng(6)
    motion_amplitude_a = rng.normal(size=(1, 6)) + 1j * rng.normal(size=(1, 6))
    motion_amplitude = np.concatenate([motion_amplitude_a, mirror * motion_amplitude_a], axis=1)
    force = near_field_mean_drift_force_with_imposed_motion(a + b, motion_amplitude, [np.pi/2]).values[0, 0, :]
    assert np.abs(force[:6]).max() > 0.0
    np.testing.assert_allclose(force[6:], mirror * force[:6], rtol=1e-8, atol=1e-8 * np.abs(force).max())


def roll_pitch_yaw_matrix(angles):
    """Rotation matrix R = R_roll R_pitch R_yaw (same convention as WAMIT)."""
    roll, pitch, yaw = angles
    R_roll = np.array([[1.0, 0.0, 0.0],
                       [0.0, np.cos(roll), -np.sin(roll)],
                       [0.0, np.sin(roll), np.cos(roll)]])
    R_pitch = np.array([[np.cos(pitch), 0.0, np.sin(pitch)],
                        [0.0, 1.0, 0.0],
                        [-np.sin(pitch), 0.0, np.cos(pitch)]])
    R_yaw = np.array([[np.cos(yaw), -np.sin(yaw), 0.0],
                      [np.sin(yaw), np.cos(yaw), 0.0],
                      [0.0, 0.0, 1.0]])
    return R_roll @ R_pitch @ R_yaw


def exact_hydrostatic_force_on_displaced_body(mesh, motion, rotation_center, rho, g):
    """Generalized forces of the six rigid body dofs (moments around the fixed `rotation_center`)
    of the hydrostatic pressure on the part below z=0 of `mesh` displaced by the finite real `motion`."""
    R = roll_pitch_yaw_matrix(motion[3:])
    vertices = (mesh.vertices - rotation_center) @ R.T + rotation_center + motion[:3]
    # 2x2 Gauss-Legendre is exact for the quadratic integrand of the moments on planar faces.
    wet = cpt.Mesh(vertices=vertices, faces=mesh.faces).immersed_part().with_quadrature("Gauss-Legendre 2")
    pressure = -rho * g * wet.quadrature_points[0][:, :, 2]
    dofs = cpt.rigid_body_dofs(rotation_center=rotation_center)
    return np.array([
        wet.surface_integral(-pressure * np.sum(dof.evaluate_motion(wet, location="quadrature_points") * wet.faces_normals[:, None, :], axis=-1))
        for dof in dofs.values()
        ])


def test_hydrostatic_force_on_moving_body_against_displaced_mesh():
    # Hydrostatic pressure on the hull and in the strip at the waterline, without waves.
    # The mean over a period of the exact force on the displaced body is compared with the diagonal of the QTF.
    # The hydrostatic part of `near_field_mean_drift_force` is reimplemented below (pressure on the hull, strip
    # pressure and width at the waterline), since the function cannot be called without the hydrodynamics.
    rho, g = 1000.0, 9.81
    rotation_center = np.array([0.4, -0.3, -0.6])
    full_mesh = cpt.mesh_parallelepiped(size=(3.0, 2.0, 2.0), resolution=(6, 4, 4), center=(0.0, 0.0, -0.5))
    mesh = full_mesh.immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=rotation_center))
    rng = np.random.default_rng(4)
    q = rng.normal(size=(1, 6)) + 1j * rng.normal(size=(1, 6))

    u, _ = body._displacement_expansions(Expansion.first_order(q))
    # Plain arrays are order-0 constants in the operations with expansions.
    hydrostatic_pressure = -rho * g * mesh.faces_centers[:, 2]
    gradient_hydrostatic_pressure = np.array([0.0, 0.0, -rho * g])
    pressure = hydrostatic_pressure + (gradient_hydrostatic_pressure * u).sum(-1)
    dofs_on_waterline = body._dofs_on_waterline()
    vertical_motion_waterline = np.sum(q[:, :, None] * dofs_on_waterline[None, :, :, 2], axis=1)
    pressure_waterline = Expansion.first_order(-rho * g * vertical_motion_waterline / 2)  # mean over the strip
    strip_width_waterline = Expansion.first_order(-vertical_motion_waterline / np.sqrt(1 - mesh.merged().normals_edges_waterline[:, 2]**2))
    force = body.integrate_pressure_on_moving_body(pressure, Expansion.first_order(q)) + body.integrate_pressure_on_waterline(pressure_waterline * strip_width_waterline)

    def mean_second_order_force(eps):
        phases = np.linspace(0.0, 2*np.pi, 32, endpoint=False)
        forces = [exact_hydrostatic_force_on_displaced_body(full_mesh, eps * np.real(q[0] * np.exp(-1j * t)), rotation_center, rho, g) for t in phases]
        still_water = exact_hydrostatic_force_on_displaced_body(full_mesh, np.zeros(6), rotation_center, rho, g)
        return (np.mean(forces, axis=0) - still_water) / eps**2  # The error is of order eps**2

    exact = mean_second_order_force(1e-3)
    np.testing.assert_allclose(force.coefficients[2][0, 0], exact, atol=1e-4 * np.max(np.abs(exact)))
