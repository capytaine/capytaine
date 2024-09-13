import pytest

import numpy as np
import xarray as xr

import capytaine as cpt
from capytaine import __version__


@pytest.fixture
def sphere():
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body


def test_exportable_settings():
    gf = cpt.Delhommeau(tabulation_nr=10, tabulation_nz=10,
                    tabulation_grid_shape="legacy", tabulation_nb_integration_points=50)
    assert gf.exportable_settings['green_function'] == 'Delhommeau'
    assert gf.exportable_settings['tabulation_nb_integration_points'] == 50
    assert gf.exportable_settings['tabulation_grid_shape'] == "legacy"
    assert gf.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'

    gf2 = cpt.XieDelhommeau()
    assert gf2.exportable_settings['green_function'] == 'XieDelhommeau'

    engine = cpt.BasicMatrixEngine(matrix_cache_size=0)
    assert engine.exportable_settings['engine'] == 'BasicMatrixEngine'
    assert engine.exportable_settings['matrix_cache_size'] == 0
    assert engine.exportable_settings['linear_solver'] == 'lu_decomposition'

    solver = cpt.BEMSolver(green_function=gf, engine=engine)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'
    assert solver.exportable_settings['engine'] == 'BasicMatrixEngine'
    assert solver.exportable_settings['matrix_cache_size'] == 0
    assert solver.exportable_settings['linear_solver'] == 'lu_decomposition'


def test_direct_solver(sphere):
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0)
    solver = cpt.BEMSolver()
    direct_result = solver.solve(problem, method='direct')
    indirect_result = solver.solve(problem, method='indirect')
    assert direct_result.forces["Surge"] == pytest.approx(indirect_result.forces["Surge"], rel=1e-1)


@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_same_result_with_symmetries(method):
    solver = cpt.BEMSolver()
    sym_mesh = cpt.ReflectionSymmetricMesh(cpt.mesh_sphere(center=(0, 2, 0)).immersed_part(), cpt.xOz_Plane)
    sym_body = cpt.FloatingBody(mesh=sym_mesh, dofs=cpt.rigid_body_dofs())
    sym_result = solver.solve(cpt.DiffractionProblem(body=sym_body, omega=1.0), method=method)
    mesh = sym_mesh.merged()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    result = solver.solve(cpt.DiffractionProblem(body=body, omega=1.0), method=method)
    assert sym_result.forces["Surge"] == pytest.approx(result.forces["Surge"], rel=1e-10)


def test_parallelization(sphere):
    pytest.importorskip("joblib")
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'radiating_dof': list(sphere.dofs.keys()),
    })
    solver.fill_dataset(test_matrix, sphere, n_jobs=2)


def test_float32_solver(sphere):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(floating_point_precision="float32"))
    pb = cpt.RadiationProblem(body=sphere, radiating_dof="Surge", omega=1.0)
    solver.solve(pb)


def test_fill_dataset(sphere):
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'wave_direction': np.linspace(0.0, np.pi, 3),
        'radiating_dof': list(sphere.dofs.keys()),
        'rho': [1025.0],
        'water_depth': [np.inf, 10.0],
        'g': [9.81]
    })
    dataset = solver.fill_dataset(test_matrix, sphere, n_jobs=1)

    # Tests on the coordinates
    assert list(dataset.coords['influenced_dof']) == list(dataset.coords['radiating_dof']) == list(sphere.dofs.keys())
    assert dataset.body_name == sphere.name
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


def test_warning_mesh_resolution(sphere, caplog):
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=sphere, wavelength=0.1*sphere.minimal_computable_wavelength)
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "resolution " in caplog.text
