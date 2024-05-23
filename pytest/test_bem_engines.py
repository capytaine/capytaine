import pytest

import numpy as np
import capytaine as cpt


def test_cache_matrices():
    """Test how the BasicMatrixEngine caches the interaction matrices."""
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part()
    gf = cpt.Delhommeau()
    params_1 = (mesh, mesh, 0.0, np.inf, 1.0, gf)
    params_2 = (mesh, mesh, 0.0, np.inf, 2.0, gf)

    # No cache
    engine = cpt.BasicMatrixEngine(matrix_cache_size=0)
    S, K             = engine.build_matrices(*params_1)
    S_again, K_again = engine.build_matrices(*params_1)
    assert S is not S_again
    assert K is not K_again

    # Cache
    engine = cpt.BasicMatrixEngine(matrix_cache_size=1)
    S, K                     = engine.build_matrices(*params_1)
    S_again, K_again         = engine.build_matrices(*params_1)
    _, _                     = engine.build_matrices(*params_2)
    S_once_more, K_once_more = engine.build_matrices(*params_2)
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

    reference_solver = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver="gmres", matrix_cache_size=0)
    )
    reference_result = reference_solver.solve(problem,method=method)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b

    my_bem_solver = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver=my_linear_solver, matrix_cache_size=0)
    )
    assert 'my_linear_solver' in my_bem_solver.exportable_settings['linear_solver']

    result = my_bem_solver.solve(problem,method=method)
    assert np.isclose(reference_result.added_masses['Surge'], result.added_masses['Surge'])
