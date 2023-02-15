import pytest
import numpy as np
import capytaine as cpt

def test_gradient_G():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.faces_centers[[12], :], mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    _, gradG_3 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    np.testing.assert_allclose(gradG_1, gradG_2[[12], :])
    for i in range(mesh.nb_faces):
        if i != 12:
            np.testing.assert_allclose(gradG_1[:, i], gradG_3[[12], i])
    np.testing.assert_allclose(gradG_1[0, 12], gradG_3[12, 12] - 0.5*mesh.faces_normals[12, :])


def test_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=True)
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[i, :]
    K__ = np.einsum('...jk,...k->...j', gradG, mesh.faces_normals)
    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)


def test_manual_reconstruction_of_normal_velocity():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=True)
    neumann_bc = np.random.rand(mesh.nb_faces)
    sources = np.linalg.solve(K, neumann_bc)
    velocities = np.einsum('ijk,j->ik', gradG, sources)
    normal_velocities = np.einsum('...k,...k->...', velocities, mesh.faces_normals)
    np.testing.assert_allclose(normal_velocities, neumann_bc)


def test_a_posteriori_scalar_product_direct_method():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, early_dot_product=False)
    D_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    D__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)
    np.testing.assert_allclose(D_, D__)

