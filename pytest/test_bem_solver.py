import pytest

import numpy as np
import xarray as xr
from numpy import pi

try:
    import joblib
except ImportError:
    joblib = None

from capytaine import __version__
from capytaine.bem.solver import BEMSolver, Nemoh
from capytaine.green_functions.delhommeau import Delhommeau, XieDelhommeau
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere

sphere = Sphere(radius=1.0, ntheta=2, nphi=3, clip_free_surface=True)
sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")


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
    assert engine.exportable_settings['linear_solver'] == 'direct'

    solver = BEMSolver(green_function=gf, engine=engine)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'
    assert solver.exportable_settings['engine'] == 'BasicMatrixEngine'
    assert solver.exportable_settings['matrix_cache_size'] == 0
    assert solver.exportable_settings['linear_solver'] == 'direct'


def test_limit_frequencies():
    """Test if how the solver answers when asked for frequency of 0 or ∞."""
    solver = BEMSolver()

    solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-1.0))

    solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-10))


@pytest.mark.skipif(joblib is None, reason='joblib is not installed')
def test_parallelization():
    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'radiating_dof': list(sphere.dofs.keys()),
    })
    solver.fill_dataset(test_matrix, sphere, n_jobs=2)


def test_fill_dataset():
    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'wave_direction': np.linspace(0.0, pi, 3),
        'radiating_dof': list(sphere.dofs.keys()),
        'rho': [1025.0],
        'water_depth': [np.infty, 10.0],
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
    assert 'incoming_waves_convention' in dataset.attrs

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
