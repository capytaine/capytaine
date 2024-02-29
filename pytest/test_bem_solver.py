import pytest

import numpy as np
import xarray as xr
from numpy import pi

try:
    import joblib
except ImportError:
    joblib = None

import capytaine as cpt
from capytaine import __version__
from capytaine.bem.solver import BEMSolver
from capytaine.green_functions.delhommeau import Delhommeau, XieDelhommeau
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere
#-----------------------------------------------------------------------------#
# Test indirect solver
#-----------------------------------------------------------------------------#
@pytest.fixture
def sphere():
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body


def test_exportable_settings():
    gf = Delhommeau(tabulation_nb_integration_points=50)
    assert gf.exportable_settings['green_function'] == 'Delhommeau'
    assert gf.exportable_settings['tabulation_nb_integration_points'] == 50
    assert gf.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'

    gf2 = XieDelhommeau()
    assert gf2.exportable_settings['green_function'] == 'XieDelhommeau'

    engine = BasicMatrixEngine(matrix_cache_size=0)
    assert engine.exportable_settings['engine'] == 'BasicMatrixEngine'
    assert engine.exportable_settings['matrix_cache_size'] == 0
    assert engine.exportable_settings['linear_solver'] == 'lu_decomposition'

    solver = BEMSolver(green_function=gf, engine=engine)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'
    assert solver.exportable_settings['engine'] == 'BasicMatrixEngine'
    assert solver.exportable_settings['matrix_cache_size'] == 0
    assert solver.exportable_settings['linear_solver'] == 'lu_decomposition'


def test_limit_frequencies(sphere):
    """Test if how the solver answers when asked for frequency of 0 or ∞."""
    solver = BEMSolver()

    solver.solve(RadiationProblem(body=sphere, omega=0.0, water_depth=np.inf))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=0.0, water_depth=1.0))

    solver.solve(RadiationProblem(body=sphere, omega=np.inf, water_depth=np.inf))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=np.inf, water_depth=10))


def test_limit_frequencies_with_symmetries():
    mesh = cpt.mesh_parallelepiped(reflection_symmetry=True).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Surge")
    pb = cpt.RadiationProblem(body=body, omega=0.0)
    solver = cpt.BEMSolver()
    res = solver.solve(pb, keep_details=True)
    assert isinstance(res.added_mass['Surge'], float)


@pytest.mark.skipif(joblib is None, reason='joblib is not installed')
def test_parallelization(sphere):
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
    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'wave_direction': np.linspace(0.0, pi, 3),
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


# TODO: move the code below to test_io_xarray.py
    # wavenumbers = wavenumber_data_array(results)
    # assert isinstance(wavenumbers, xr.DataArray)

#-----------------------------------------------------------------------------#
# Test direct solver
#-----------------------------------------------------------------------------#
def test_direct_solver(sphere):
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0)
    solver = cpt.BEMSolver()
    direct_result = solver.solve(problem, method='direct')
    indirect_result = solver.solve(problem, method='indirect')
    assert direct_result.forces["Surge"] == pytest.approx(indirect_result.forces["Surge"], rel=1e-1)
