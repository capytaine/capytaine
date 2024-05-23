"""Tests for the cpt.Delhommeau() class"""

import pytest

import numpy as np
import capytaine as cpt


def test_no_tabulation():
    mesh = cpt.mesh_sphere().immersed_part()
    tabed_gf = cpt.Delhommeau()
    untabed_gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0)
    assert np.allclose(untabed_gf.evaluate(mesh, mesh)[0], tabed_gf.evaluate(mesh, mesh)[0], atol=1e-2)


# def test_panels_near_free_surface():
#     area = 4.0
#     def gf_at_depth(d):
#         mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), resolution=(1, 1), center=(0, 0, d))
#         x = np.array([[0, 0, 0]])
#         S, _ = cpt.Delhommeau().evaluate(x, mesh, free_surface=0, water_depth=np.infty, wavenumber=3.0)
#         return S
#     gf_at_depth(0.0)
#     gf_at_depth(1e-9)
#     gf_at_depth(1e-5)
#     gf_at_depth(1e-2)


def test_exact_integration_of_rankine_terms():
    gf = cpt.Delhommeau()
    center = np.array([0.0, 0.0, 0.0])
    area = 1.0

    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(1, 1))
    rankine_g_once = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(11, 11))
    # Use odd number of panels to avoid having the center on a corner of a panel, which is not defined for the strong singularity of the derivative.
    rankine_g_parts = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    assert rankine_g_once == pytest.approx(rankine_g_parts.sum(), abs=1e-2)


def test_exact_integration_with_nb_integration_points():
    # Test that the tabulation_nb_integration_points is actually taken into account.
    mesh = cpt.mesh_rectangle(center=(0, 0, -1), resolution=(1, 1))
    point = np.array([[100, 0, -1]])
    gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nb_integration_points=11)
    val1 = gf.evaluate(point, mesh)
    gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nb_integration_points=101)
    val2 = gf.evaluate(point, mesh)
    assert np.abs(val1[0] - val2[0]) > 1e-1


def test_low_freq_with_rankine_part_singularities():
    import numpy as np
    import capytaine as cpt
    gf1 = cpt.Delhommeau(gf_singularities="low_freq")
    gf2 = cpt.Delhommeau(gf_singularities="low_freq_with_rankine_part")
    point = np.array([[0.0, 0.0, -1.0]])
    area = 1.0
    wavenumber = 1.0
    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=(0, 0, -1), resolution=(1, 1))
    S1, K1 = gf1.evaluate(point, mesh, wavenumber=wavenumber, early_dot_product=False)
    S2, K2 = gf2.evaluate(point, mesh, wavenumber=wavenumber, early_dot_product=False)
    assert S1 == S2
    assert np.allclose(K1[:, :, :2], K2[:, :, :2], atol=1e-12)
    assert np.allclose(K1[:, :, 2].imag, K2[:, :, 2].imag, atol=1e-12)
    assert not np.allclose(K1[0, 0, 2].real, K2[0, 0, 2].real, atol=1e-12)
    assert np.allclose(K1[0, 0, 2].real, K2[0, 0, 2].real, atol=1e-1)


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
