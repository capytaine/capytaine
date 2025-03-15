import pytest

import numpy as np
import capytaine as cpt

from capytaine.tools.prony_decomposition import exponential_decomposition, find_best_exponential_decomposition

list_of_faces = [
    cpt.Mesh(vertices=[[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [1.0, 1.0, -1.0], [0.0, 1.0, -1.0]], faces=np.array([[0, 1, 2, 3]])),
    cpt.Mesh(vertices=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], faces=np.array([[0, 1, 2, 3]])),
        ]

@pytest.mark.parametrize("face", list_of_faces, ids=["immersed", "free_surface"])
def test_deep_water_asymptotics(face):
    gf = cpt.Delhommeau()
    k = 1.0
    depth = 1000.0
    ambda, a, nexp = gf.fortran_core.old_prony_decomposition.lisc(k*depth*np.tanh(k*depth), k*depth)
    s_inf, k_inf = gf.fortran_core.green_wave.integral_of_wave_part_infinite_depth(
        face.faces_centers[0, :], face.faces_centers[0, :], face.faces_areas[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, *gf.all_tabulation_parameters, gf.gf_singularities_index, True
    )
    s_finite, k_finite = gf.fortran_core.green_wave.integral_of_wave_part_finite_depth(
        face.faces_centers[0, :], face.vertices[face.faces[0, :], :], face.faces_centers[0, :] , face.faces_normals[0, :],
        face.faces_areas[0], face.faces_radiuses[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, depth, *gf.all_tabulation_parameters, ambda, a, True
    )
    np.testing.assert_allclose(s_inf, s_finite, rtol=1e-2)
    np.testing.assert_allclose(k_inf, k_finite, rtol=1e-2)


def test_prony_decomposition():
    x = np.linspace(0.0, 1.0, 100)
    y = 2*np.exp(-x) - 4*np.exp(-2*x)
    a, lamda = exponential_decomposition(x, y, 2)
    assert ((np.allclose(a, [-4.0, 2.0]) and np.allclose(lamda, [-2.0, -1.0]))
            or (np.allclose(a, [2.0, -4.0]) and np.allclose(lamda, [-1.0, -2.0])))


def test_find_best_prony_decomposition():
    rng = np.random.default_rng(seed=987654321)
    # May fail for other RNGs, for instance if two values from lamda_refs are very close from each other
    n_ref = 6
    a_ref = rng.uniform(-10.0, 10.0, size=n_ref)
    lamda_ref = -rng.chisquare(1.0, size=n_ref)
    def f(x):
        return sum(a_ref[i]*np.exp(lamda_ref[i]*x) for i in range(n_ref))
    a, lamda = find_best_exponential_decomposition(f, 0.0, 100.0, range(2, 15), tol=1e-4)
    assert np.allclose(np.sort(a), np.sort(a_ref))
    assert np.allclose(np.sort(lamda_ref), np.sort(lamda_ref))


def test_python_and_fortran_prony_decomposition_for_green_function():
    gf = cpt.Delhommeau(finite_depth_prony_decomposition_method="python")
    decomp_default = gf.find_best_exponential_decomposition(1.0)
    decomp_f = gf.find_best_exponential_decomposition(1.0, method="fortran")
    decomp_p = gf.find_best_exponential_decomposition(1.0, method="python")
    assert np.allclose(decomp_default[0], decomp_p[0])
    assert np.allclose(decomp_p[0], decomp_f[0], rtol=0.2)
    assert np.allclose(decomp_p[1], decomp_f[1], rtol=0.2)


def test_fingreen3D():
    k = 5.0
    h = 2.0
    mesh = cpt.mesh_sphere(radius=0.5, center=(0.0, 0.0, -1.0), resolution=(4, 4))
    S, D = cpt.FinGreen3D().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    S_ref, D_ref = cpt.Delhommeau().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    np.testing.assert_allclose(S, S_ref, rtol=1e-1, atol=1e-1)
    np.testing.assert_allclose(D, D_ref, rtol=1e-1, atol=1e-1)

@pytest.mark.parametrize("depth", [0.5, 1.0, 2.0])
def test_fingreen3D_roots(depth):
    from scipy.optimize import newton
    nk = 100
    omega2_over_g = 1.0
    all_roots = cpt.FinGreen3D().fortran_core.green_wave.dispersion_roots(nk, omega2_over_g, depth)
    wavenumber = newton(lambda x: omega2_over_g - x*np.tanh(x*depth), omega2_over_g)
    assert np.isclose(wavenumber, all_roots[0], rtol=1e-3)
    roots_idempotent = newton(lambda x: -omega2_over_g - x*np.tan(x*depth), all_roots[1:])
    np.testing.assert_allclose(roots_idempotent, all_roots[1:], rtol=1e-3)

    roots_ = newton(lambda x: -omega2_over_g - x*np.tan(x*depth), (np.arange(1, nk) + 0.25)*np.pi/depth)
    np.testing.assert_allclose(roots_, all_roots[1:], rtol=1e-4)

