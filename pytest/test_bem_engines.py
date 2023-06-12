import pytest

import numpy as np
import capytaine as cpt

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
    params_1 = (sphere.mesh, sphere.mesh, 0.0, np.infty, 1.0, gf)
    params_2 = (sphere.mesh, sphere.mesh, 0.0, np.infty, 2.0, gf)

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
    problem = RadiationProblem(body=sphere, omega=1.0, water_depth=np.infty)

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


def test_gradiant_G_shape():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    S, gradG_1 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    assert gradG_1.shape == (mesh.nb_faces, mesh.nb_faces, 3)


def test_gradient_G_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=True)

    # Scalar product with the faces normal vectors
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[i, :]

    # Shorter implementiation of the exact same scalar products
    K__ = np.einsum('...jk,...k->...j', gradG, mesh.faces_normals)

    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)


def test_gradient_G_with_collocation_points():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.faces_centers, mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    np.testing.assert_allclose(gradG_1, gradG_2)


def test_gradient_G_diagonal_term():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    # Passing two different Python objects
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    # Passing the same Python object
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False)

    diag_normal = np.zeros_like(gradG_2)
    for i in range(mesh.nb_faces):
        diag_normal[i, i, :] = mesh.faces_normals[i, :]

    np.testing.assert_allclose(gradG_1, gradG_2 - 0.5*diag_normal)


def test_a_posteriori_scalar_product_direct_method():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False)
    D_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    D__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)
    np.testing.assert_allclose(D_, D__)

