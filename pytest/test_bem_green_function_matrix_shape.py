"""Tests for the matrices returned by cpt.Delhommeau() class"""

import pytest

import numpy as np
import capytaine as cpt


@pytest.mark.parametrize("adjoint_double_layer", [True, False])
def test_gradiant_G_shape(adjoint_double_layer):
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    S, gradG_1 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=False,adjoint_double_layer=adjoint_double_layer)
    assert gradG_1.shape == (mesh.nb_faces, mesh.nb_faces, 3)


def test_gradient_G_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=True)

    # Scalar product with the faces normal vectors
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[i, :]

    # Shorter implementiation of the exact same scalar products
    K__ = np.einsum('...jk,...k->...j', gradG, mesh.faces_normals)

    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)


def test_gradient_G_a_posteriori_scalar_product_directBIE():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=False, adjoint_double_layer=False)
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=True, adjoint_double_layer=False)

    # Scalar product with the faces normal vectors
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]

    # Shorter implementiation of the exact same scalar products
    K__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)

    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)
    np.testing.assert_allclose(K_, K__)


def test_gradient_G_with_collocation_points():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.faces_centers, mesh, 0.0, np.inf, 1.0, early_dot_product=False)
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, np.inf, 1.0, early_dot_product=False)
    np.testing.assert_allclose(gradG_1, gradG_2, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("adjoint_double_layer", [True, False])
def test_gradient_G_diagonal_term(adjoint_double_layer):
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    # Passing two different Python objects
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, np.inf, 1.0,
                                           early_dot_product=False, adjoint_double_layer=adjoint_double_layer)
    # Passing the same Python object
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0,
                                           early_dot_product=False, adjoint_double_layer=adjoint_double_layer)

    diag_normal = np.zeros_like(gradG_2)
    for i in range(mesh.nb_faces):
        diag_normal[i, i, :] = mesh.faces_normals[i, :]

    np.testing.assert_allclose(gradG_1, gradG_2 - 0.5*diag_normal, atol=1e-10, rtol=1e-10)


def test_a_posteriori_scalar_product_direct_method():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.inf, 1.0, early_dot_product=False)
    D_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    D__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)
    np.testing.assert_allclose(D_, D__)
