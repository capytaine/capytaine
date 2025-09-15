import pytest

import numpy as np
import capytaine as cpt


def test_engine_repr():
    r = repr(cpt.BasicMatrixEngine())
    assert "Delhommeau" in r


def test_cache_matrices():
    """Test how the BasicMatrixEngine caches the interaction matrices."""
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part()
    params_1 = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, adjoint_double_layer=True)
    params_2 = dict(free_surface=0.0, water_depth=np.inf, wavenumber=2.0, adjoint_double_layer=True)

    # No cache
    engine = cpt.BasicMatrixEngine()
    S, K             = engine.build_matrices_with_symmetries(mesh, mesh, **params_1)
    S_again, K_again = engine.build_matrices_with_symmetries(mesh, mesh, **params_1)
    assert S is not S_again
    assert K is not K_again

    # Cache
    engine = cpt.BasicMatrixEngine()
    S, K                     = engine.build_and_cache_matrices_with_symmetries(mesh, mesh, **params_1)
    S_again, K_again         = engine.build_and_cache_matrices_with_symmetries(mesh, mesh, **params_1)
    _, _                     = engine.build_and_cache_matrices_with_symmetries(mesh, mesh, **params_2)
    S_once_more, K_once_more = engine.build_and_cache_matrices_with_symmetries(mesh, mesh, **params_2)
    assert S is S_again
    assert S is not S_once_more
    assert K is K_again
    assert K is not K_once_more


@pytest.mark.parametrize("method", ['indirect', 'direct'])
def test_custom_linear_solver(method):
    """Solve a simple problem with a custom linear solver."""
    sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part())
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)

    reference_solver = cpt.BEMSolver(engine=cpt.BasicMatrixEngine(), method=method)
    reference_result = reference_solver.solve(problem)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b

    my_bem_solver = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver=my_linear_solver),
        method=method
    )
    assert 'my_linear_solver' in my_bem_solver.exportable_settings['linear_solver']

    result = my_bem_solver.solve(problem)
    assert np.isclose(reference_result.added_masses['Surge'], result.added_masses['Surge'])


@pytest.mark.parametrize("method", ['indirect', 'direct'])
def test_custom_linear_solver_returning_wrong_shape(method):
    sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part())
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b.reshape(-1, 1)

    my_bem_solver = cpt.BEMSolver(
        method=method,
        engine=cpt.BasicMatrixEngine(linear_solver=my_linear_solver)
    )
    with pytest.raises(ValueError):
        my_bem_solver.solve(problem)
