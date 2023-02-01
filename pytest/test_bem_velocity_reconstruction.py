import pytest
import numpy as np
import capytaine as cpt

def test_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, K_dim=3)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, K_dim=1)
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[i, :]
    K__ = np.einsum('...jk,...k->...j', gradG, mesh.faces_normals)
    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)


def test_a_posteriori_scalar_product_direct_method():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, -np.infty, 1.0, K_dim=3)
    D_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    D__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)
    np.testing.assert_allclose(D_, D__)
