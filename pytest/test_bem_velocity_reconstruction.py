import pytest
import numpy as np
import capytaine as cpt


def test_gradiant_G_shape():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    S, gradG_1 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    assert gradG_1.shape == (mesh.nb_faces, mesh.nb_faces, 3)


def test_gradient_G_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=True)

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
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.faces_centers, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    np.testing.assert_allclose(gradG_1, gradG_2)


def test_gradient_G_diagonal_term():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    # Passing two different Python objects
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    # Passing the same Python object
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)

    diag_normal = np.zeros_like(gradG_2)
    for i in range(mesh.nb_faces):
        diag_normal[i, i, :] = mesh.faces_normals[i, :]

    np.testing.assert_allclose(gradG_1, gradG_2 - 0.5*diag_normal)


@pytest.fixture
def solver():
    return cpt.BEMSolver()

@pytest.fixture
def solved_problem(solver):
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs={"random_dof": np.random.rand(mesh.nb_faces, 3)})
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="random_dof")
    res = solver.solve(pb)
    return res


def test_velocity_at_point(solver, solved_problem):
    point = np.array([[10.0, 10.0, -2.0]])
    velocities = solver.get_velocity(solved_problem, point)
    assert velocities.shape == (1, 3)


def test_reconstruction_of_given_boundary_condition(solver, solved_problem):
    velocities = solver.get_velocity(solved_problem, solved_problem.body.mesh)
    normal_velocities = np.einsum('...k,...k->...', velocities, solved_problem.body.mesh.faces_normals)
    np.testing.assert_allclose(normal_velocities, solved_problem.problem.boundary_condition)


def test_a_posteriori_scalar_product_direct_method():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    D_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    D__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)
    np.testing.assert_allclose(D_, D__)

