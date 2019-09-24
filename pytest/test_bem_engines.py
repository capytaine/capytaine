import pytest

import numpy as np

from capytaine.bem.solver import BEMSolver
from capytaine.green_functions.delhommeau import Delhommeau
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere

sphere = Sphere(radius=1.0, ntheta=2, nphi=3, clip_free_surface=True)
sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")

def test_cache_matrices():
    """Test how the BasicMatrixEngine caches the interaction matrices."""
    gf = Delhommeau()
    params_1 = (sphere.mesh, sphere.mesh, 0.0, -np.infty, 1.0, gf)
    params_2 = (sphere.mesh, sphere.mesh, 0.0, -np.infty, 2.0, gf)

    # No cache
    engine = BasicMatrixEngine(matrix_cache_size=0)
    S, K             = engine.build_matrices(*params_1)
    S_again, K_again = engine.build_matrices(*params_1)
    assert S is not S_again
    assert K is not K_again

    # Cache
    engine = BasicMatrixEngine(matrix_cache_size=1)
    S, K                     = engine.build_matrices(*params_1)
    S_again, K_again         = engine.build_matrices(*params_1)
    _, _                     = engine.build_matrices(*params_2)
    S_once_more, K_once_more = engine.build_matrices(*params_2)
    assert S is S_again
    assert S is not S_once_more
    assert K is K_again
    assert K is not K_once_more


def test_custom_linear_solver():
    """Solve a simple problem with a custom linear solver."""
    problem = RadiationProblem(body=sphere, omega=1.0, sea_bottom=-np.infty)

    reference_solver = BEMSolver(
        engine=BasicMatrixEngine(linear_solver="gmres", matrix_cache_size=0)
    )
    reference_result = reference_solver.solve(problem)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b

    my_bem_solver = BEMSolver(
        engine=BasicMatrixEngine(linear_solver=my_linear_solver, matrix_cache_size=0)
    )
    assert 'my_linear_solver' in my_bem_solver.exportable_settings['linear_solver']

    result = my_bem_solver.solve(problem)
    assert np.isclose(reference_result.added_masses['Surge'], result.added_masses['Surge'])

