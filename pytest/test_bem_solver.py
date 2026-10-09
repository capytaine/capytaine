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
import pytest

import numpy as np
import xarray as xr

import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine import __version__

from capytaine.meshes.symmetric_meshes import ReflectionSymmetricMesh

def test_exportable_settings():
    gf = cpt.Delhommeau(
            tabulation_nr=10, tabulation_nz=10,
            tabulation_grid_shape="legacy",
            tabulation_nb_integration_points=50,
            finite_depth_prony_decomposition_method="fortran"
            )
    assert gf.exportable_settings['green_function'] == 'Delhommeau'
    assert gf.exportable_settings['tabulation_nb_integration_points'] == 50
    assert gf.exportable_settings['tabulation_grid_shape'] == "legacy"
    assert gf.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'

    gf2 = cpt.XieDelhommeau()
    assert gf2.exportable_settings['green_function'] == 'XieDelhommeau'

    engine = cpt.DefaultMatrixEngine(green_function=gf)
    assert engine.exportable_settings['engine'] == 'DefaultMatrixEngine'
    assert engine.exportable_settings['linear_solver'] == 'lu_decomposition'

    solver = cpt.BEMSolver(engine=engine)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'
    assert solver.exportable_settings['engine'] == 'DefaultMatrixEngine'
    assert solver.exportable_settings['linear_solver'] == 'lu_decomposition'

    solver = cpt.BEMSolver(green_function=gf)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'

def test_cannot_define_gf_and_engine_in_solver():
    with pytest.raises(ValueError):
        cpt.BEMSolver(engine=cpt.DefaultMatrixEngine(), green_function=cpt.Delhommeau())

def test_solver_has_initialized_timer():
    s = helpers.solver()
    assert s.timer.total == 0.0

def test_solver_update_timer():
    sphere = helpers.small_sphere_body()
    solver = helpers.solver()
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0)
    solver.solve(problem)
    assert solver.timer.total > 0.0

def test_direct_solver():
    sphere = helpers.small_sphere_body()
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0)
    direct_solver = helpers.solver(method='direct')
    direct_result = direct_solver.solve(problem)
    indirect_solver = helpers.solver(method='indirect')
    indirect_result = indirect_solver.solve(problem)
    assert direct_result.forces["Surge"] == pytest.approx(indirect_result.forces["Surge"], rel=1e-1)


@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_same_result_with_symmetries(method):
    solver = helpers.solver(method=method)
    sym_mesh = ReflectionSymmetricMesh(cpt.mesh_sphere(center=(0, 2, 0)).immersed_part(), plane='xOz')
    sym_body = cpt.FloatingBody(mesh=sym_mesh, dofs=cpt.rigid_body_dofs())
    sym_result = solver.solve(cpt.DiffractionProblem(body=sym_body, omega=1.0))
    mesh = sym_mesh.merged()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    result = solver.solve(cpt.DiffractionProblem(body=body, omega=1.0))
    assert sym_result.forces["Surge"] == pytest.approx(result.forces["Surge"], rel=1e-4)


def test_parallelization():
    sphere = helpers.small_sphere_body()
    pytest.importorskip("joblib")
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'radiating_dof': list(sphere.dofs.keys()),
    })
    solver.fill_dataset(test_matrix, sphere, n_jobs=2)


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("n_threads", [1, 2])
def test_control_threads(n_jobs, n_threads):
    sphere = helpers.small_sphere_body()
    pytest.importorskip("joblib")
    pytest.importorskip("threadpoolctl")
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'radiating_dof': list(sphere.dofs.keys()),
    })
    solver.fill_dataset(test_matrix, sphere, n_jobs=n_jobs, n_threads=n_threads)


def test_nb_timer():
    sphere = helpers.small_sphere_body()
    pytest.importorskip("joblib")
    from joblib import cpu_count
    solver = helpers.solver()
    n_jobs = min(cpu_count(), 3)
    problems = [
            cpt.RadiationProblem(body=sphere, radiating_dof="Surge", omega=omega)
            for omega in np.linspace(0.1, 3.0, 5)
            ]
    solver.solve_all(problems, n_jobs=n_jobs)
    assert len(solver.timer_summary().columns) == n_jobs


def test_float32_solver():
    sphere = helpers.small_sphere_body()
    solver = helpers.solver(floating_point_precision="float32")
    pb = cpt.RadiationProblem(body=sphere, radiating_dof="Surge", omega=1.0)
    result = solver.solve(pb)
    assert result.pressure.dtype == 'complex64' and result.potential.dtype == 'complex64'


def test_LiangWuNoblesseGF():
    sphere = helpers.small_sphere_body()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'radiating_dof': list(sphere.dofs),
    })
    solver = cpt.BEMSolver(green_function=cpt.LiangWuNoblesseGF())
    ref_solver = helpers.solver()
    ds = solver.fill_dataset(test_matrix, sphere)
    ref_ds = ref_solver.fill_dataset(test_matrix, sphere)
    assert np.allclose(ds.added_mass.values, ref_ds.added_mass.values, rtol=1e-2)


def test_fill_dataset():
    sphere = helpers.small_sphere_body()
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'wave_direction': np.linspace(0.0, np.pi, 3),
        'radiating_dof': list(sphere.dofs.keys()),
        'rho': [1025.0],
        'water_depth': [np.inf, 30.0],
        'g': [9.81]
    })
    dataset = solver.fill_dataset(test_matrix, sphere, n_jobs=1)

    # Tests on the coordinates
    assert list(dataset.coords['influenced_dof']) == list(dataset.coords['radiating_dof']) == list(sphere.dofs.keys())
    assert dataset.rho == test_matrix.rho
    assert dataset.g == test_matrix.g

    # Tests on the results
    assert 'added_mass' in dataset
    assert 'radiation_damping' in dataset
    assert 'Froude_Krylov_force' in dataset
    assert 'diffraction_force' in dataset

    # Test the attributes
    assert dataset.attrs['capytaine_version'] == __version__
    assert 'start_of_computation' in dataset.attrs

    # Try to strip out the outputs and recompute
    naked_data = dataset.drop_vars(["added_mass", "radiation_damping", "diffraction_force", "Froude_Krylov_force"])
    recomputed_dataset = solver.fill_dataset(naked_data, [sphere])
    assert recomputed_dataset.rho == dataset.rho
    assert recomputed_dataset.g == dataset.g
    assert "added_mass" in recomputed_dataset
    assert np.allclose(recomputed_dataset["added_mass"].data, dataset["added_mass"].data)


def test_warning_mesh_resolution(caplog):
    sphere = helpers.small_sphere_body()
    solver = helpers.solver()
    pb = cpt.RadiationProblem(body=sphere, wavelength=0.1*sphere.minimal_computable_wavelength)
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "resolution " in caplog.text
