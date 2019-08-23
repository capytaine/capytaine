import pytest

import numpy as np
import xarray as xr
from numpy import pi

from capytaine import __version__
from capytaine.bem.solver import BEMSolver, Nemoh
from capytaine.bem.green_functions import Delhommeau
from capytaine.bem.engines import BasicEngine
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere

sphere = Sphere(radius=1.0, ntheta=2, nphi=3, clip_free_surface=True)
sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")


def test_exportable_settings():
    gf = Delhommeau(tabulation_nb_integration_points=50)
    assert gf.exportable_settings['green_function'] == 'Delhommeau'
    assert gf.exportable_settings['tabulation_nb_integration_points'] == 50
    assert gf.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'

    engine = BasicEngine(matrix_cache_size=0)
    assert engine.exportable_settings['engine'] == 'BasicEngine'
    assert engine.exportable_settings['matrix_cache_size'] == 0
    assert engine.exportable_settings['linear_solver'] == 'gmres'

    solver = BEMSolver(green_functions=gf, engine=engine)
    assert solver.exportable_settings['green_function'] == 'Delhommeau'
    assert solver.exportable_settings['tabulation_nb_integration_points'] == 50
    assert solver.exportable_settings['finite_depth_prony_decomposition_method'] == 'fortran'
    assert solver.exportable_settings['engine'] == 'BasicEngine'
    assert solver.exportable_settings['matrix_cache_size'] == 0
    assert solver.exportable_settings['linear_solver'] == 'gmres'


def test_cache_matrices():
    """Test how the solver caches the interaction matrices."""
    params_1 = dict(free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0)
    params_2 = dict(free_surface=0.0, sea_bottom=-np.infty, wavenumber=2.0)

    # No cache
    solver = Nemoh(matrix_cache_size=0)
    S, K = solver.build_matrices(sphere.mesh, sphere.mesh, **params_1)
    S_again, K_again = solver.build_matrices(sphere.mesh, sphere.mesh, **params_1)
    assert S is not S_again
    assert K is not K_again

    # Cache
    solver = Nemoh(matrix_cache_size=1)
    S, K = solver.build_matrices(sphere.mesh, sphere.mesh, **params_1)
    S_again, K_again = solver.build_matrices(sphere.mesh, sphere.mesh, **params_1)
    _, _ = solver.build_matrices(sphere.mesh, sphere.mesh, **params_2)
    S_once_more, K_once_more = solver.build_matrices(sphere.mesh, sphere.mesh, **params_1)
    assert S is S_again
    assert S is not S_once_more
    assert K is K_again
    assert K is not K_once_more


def test_cache_rankine_matrices():
    """Test how the solver caches the Rankine part of the interaction matrices."""
    params_1 = dict(free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0)

    # No cache of rankine matrices
    solver = Nemoh(matrix_cache_size=1, cache_rankine_matrices=True)
    Sr, Vr = solver.build_matrices_rankine(sphere.mesh, sphere.mesh, **params_1)
    Sr_again, Vr_again = solver.build_matrices_rankine(sphere.mesh, sphere.mesh, **params_1)
    assert Sr is Sr_again
    assert Vr is Vr_again

    # Cache of rankine matrices
    solver = Nemoh(matrix_cache_size=1)
    Sr, Vr = solver.build_matrices_rankine(sphere.mesh, sphere.mesh, **params_1)
    Sr_again, Vr_again = solver.build_matrices_rankine(sphere.mesh, sphere.mesh, **params_1)
    assert Sr is not Sr_again
    assert Vr is not Vr_again


def test_custom_linear_solver():
    """Solve a simple problem with a custom linear solver."""
    problem = RadiationProblem(body=sphere, omega=1.0, sea_bottom=-np.infty)

    reference_solver = Nemoh(linear_solver="gmres", matrix_cache_size=0, hierarchical_matrices=False)
    reference_result = reference_solver.solve(problem)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b

    my_bem_solver = Nemoh(linear_solver=my_linear_solver, matrix_cache_size=0, hierarchical_matrices=False)
    assert 'my_linear_solver' in my_bem_solver.exportable_settings()['linear_solver']

    result = my_bem_solver.solve(problem)
    assert np.isclose(reference_result.added_masses['Surge'], result.added_masses['Surge'])


def test_limit_frequencies():
    """Test if how the solver answers when asked for frequency of 0 or ∞."""
    solver = Nemoh()

    solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-1.0))

    solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-10))


def test_fill_dataset():
    solver = Nemoh()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4.0, 3),
        'wave_direction': np.linspace(0.0, pi, 3),
        'radiating_dof': list(sphere.dofs.keys()),
        'rho': [1025.0],
        'water_depth': [np.infty, 10.0]
    })
    dataset = solver.fill_dataset(test_matrix, [sphere])

    # Tests on the coordinates
    assert list(dataset.coords['influenced_dof']) == list(dataset.coords['radiating_dof']) == list(sphere.dofs.keys())
    assert dataset.body_name == sphere.name
    assert dataset.rho == test_matrix.rho

    # Tests on the results
    assert 'added_mass' in dataset
    assert 'radiation_damping' in dataset
    assert 'Froude_Krylov_force' in dataset
    assert 'diffraction_force' in dataset

    # Test the attributes
    assert dataset.attrs['capytaine_version'] == __version__
    assert 'start_of_computation' in dataset.attrs
    assert 'cache_rankine_matrices' in dataset.attrs
    assert 'incoming_waves_convention' in dataset.attrs

    # Try to strip out the outputs and recompute
    naked_data = dataset.drop(["added_mass", "radiation_damping", "diffraction_force", "Froude_Krylov_force"])
    recomputed_dataset = solver.fill_dataset(naked_data, [sphere])
    assert recomputed_dataset.rho == dataset.rho
    assert "added_mass" in recomputed_dataset
    assert np.allclose(recomputed_dataset["added_mass"].data, dataset["added_mass"].data)


def test_fill_dataset_with_kochin_functions():
    solver = Nemoh()
    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': np.linspace(0, 2*pi, 5),
        'radiating_dof': list(sphere.dofs.keys()),
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin' in ds


# TODO: move the code below to test_io_xarray.py
    # wavenumbers = wavenumber_data_array(results)
    # assert isinstance(wavenumbers, xr.DataArray)
