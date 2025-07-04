import pytest

import numpy as np
import capytaine as cpt

from capytaine.tools.prony_decomposition import exponential_decomposition, find_best_exponential_decomposition

faces_examples = {
        "immersed": cpt.Mesh(vertices=[[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [1.0, 1.0, -1.0], [0.0, 1.0, -1.0]], faces=np.array([[0, 1, 2, 3]])),
        "free_surface": cpt.Mesh(vertices=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]], faces=np.array([[0, 1, 2, 3]])),
        }

@pytest.mark.parametrize("derivative_with_respect_to_first_variable", [True, False])
@pytest.mark.parametrize("face_location", ["immersed", "free_surface"])
def test_deep_water_asymptotics_of_wave_terms_low_freq(derivative_with_respect_to_first_variable, face_location):
    gf = cpt.Delhommeau()
    gf_singularities_index = gf.gf_singularities_fortran_enum["low_freq"]
    face = faces_examples[face_location]
    k = 1.0
    depth = 1000.0
    s_inf, k_inf = gf.fortran_core.green_wave.integral_of_wave_part_infinite_depth(
        face.faces_centers[0, :], face.faces_centers[0, :], face.faces_areas[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, *gf.all_tabulation_parameters, gf_singularities_index, derivative_with_respect_to_first_variable
    )
    s_finite, k_finite = gf.fortran_core.green_wave.integral_of_wave_parts_finite_depth(
        face.faces_centers[0, :], face.faces_centers[0, :], face.faces_areas[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, depth, *gf.all_tabulation_parameters, gf_singularities_index, derivative_with_respect_to_first_variable
    )
    np.testing.assert_allclose(s_inf, s_finite, rtol=1e-2)
    np.testing.assert_allclose(k_inf, k_finite, rtol=1e-2)


@pytest.mark.parametrize("derivative_with_respect_to_first_variable", [True, False])
@pytest.mark.parametrize("face_location", ["immersed"])
def test_deep_water_asymptotics_of_wave_terms_high_freq(derivative_with_respect_to_first_variable, face_location):
    gf = cpt.Delhommeau()
    gf_singularities_index = gf.gf_singularities_fortran_enum["high_freq"]
    face = faces_examples[face_location]
    k = 1.0
    depth = 1000.0
    s_inf, k_inf = gf.fortran_core.green_wave.integral_of_wave_part_infinite_depth(
        face.faces_centers[0, :], face.faces_centers[0, :], face.faces_areas[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, *gf.all_tabulation_parameters, gf_singularities_index, derivative_with_respect_to_first_variable
    )
    s_finite, k_finite = gf.fortran_core.green_wave.integral_of_wave_parts_finite_depth(
        face.faces_centers[0, :], face.faces_centers[0, :], face.faces_areas[0], face.quadrature_points[0][0, :, :], face.quadrature_points[1][0],
        k, depth, *gf.all_tabulation_parameters, gf_singularities_index, derivative_with_respect_to_first_variable
    )
    np.testing.assert_allclose(s_inf, s_finite, rtol=1e-2)
    np.testing.assert_allclose(k_inf, k_finite, rtol=1e-2)


@pytest.mark.parametrize("derivative_with_respect_to_first_variable", [True, False])
@pytest.mark.parametrize("face_location", ["immersed", "free_surface"])
def test_deep_water_asymptotics_of_prony_decomposition(derivative_with_respect_to_first_variable, face_location):
    gf = cpt.Delhommeau()
    k = 1.0
    depth = 1000.0
    face = faces_examples[face_location]
    decomp = gf.find_best_exponential_decomposition(k*depth, method="python")
    s_prony, k_prony = gf.fortran_core.green_wave.integral_of_prony_decomp_finite_depth(
            face.faces_centers[0, :], face.vertices[face.faces[0, :], :], face.faces_centers[0, :], face.faces_normals[0, :], face.faces_areas[0], face.faces_radiuses[0],
            depth, decomp, derivative_with_respect_to_first_variable)
    np.testing.assert_allclose(s_prony, 0.0, atol=1e-3)
    np.testing.assert_allclose(k_prony, 0.0, atol=1e-3)


def test_high_freq_vs_low_freq():
    gf = cpt.Delhommeau()
    mesh = cpt.mesh_sphere().immersed_part()
    k = 1.0
    depth = 10.0
    decomp = gf.find_best_exponential_decomposition(k*depth, method="python")
    def compute_matrices(gf_singularities):
        adjoint_double_layer = True
        S, K = gf._init_matrices(
            (mesh.nb_faces, mesh.nb_faces), early_dot_product=True
        )
        gf.fortran_core.matrices.build_matrices(
                mesh.faces_centers, mesh.faces_normals,
                mesh.vertices, mesh.faces + 1,
                mesh.faces_centers, mesh.faces_normals,
                mesh.faces_areas, mesh.faces_radiuses,
                *mesh.quadrature_points,
                k, depth,
                *gf.all_tabulation_parameters,
                gf.finite_depth_method_index, decomp,
                gf.dispersion_relation_roots,
                gf.gf_singularities_fortran_enum[gf_singularities],
                adjoint_double_layer,
                S, K,
                )
        return S, K
    S_low, K_low = compute_matrices('low_freq')
    S_high, K_high = compute_matrices('high_freq')
    np.testing.assert_allclose(S_low, S_high, atol=1e-2*np.linalg.norm(S_low))
    np.testing.assert_allclose(K_low, K_high, atol=2e-2*np.linalg.norm(K_low))



def test_infinite_frequency_prony_decomposition():
    gf = cpt.Delhommeau()
    with pytest.raises(NotImplementedError):
        gf.find_best_exponential_decomposition(np.inf, method='fortran')
    gf.find_best_exponential_decomposition(np.inf, method='python')


def test_infinite_frequency():
    gf = cpt.Delhommeau(gf_singularities='high_freq')
    mesh = cpt.mesh_sphere().immersed_part()
    S_, K_ = gf.evaluate(mesh, mesh, free_surface=0.0, water_depth=10.0, wavenumber=1000.0)
    S, K = gf.evaluate(mesh, mesh, free_surface=0.0, water_depth=10.0, wavenumber=np.inf)
    np.testing.assert_allclose(S/np.linalg.norm(S), S_/np.linalg.norm(S_), atol=1e-3)
    np.testing.assert_allclose(K/np.linalg.norm(K), K_/np.linalg.norm(K_), atol=1e-3)


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
    a, lamda = find_best_exponential_decomposition(f, 0.0, 100.0, range(2, 15), tol=1e-4, noise_on_domain_points_std=0.0)
    assert np.allclose(np.sort(a), np.sort(a_ref))
    assert np.allclose(np.sort(lamda_ref), np.sort(lamda_ref))


def test_failure_for_low_kh():
    gf = cpt.Delhommeau()

    # Legacy method returning crappy value
    gf.find_best_exponential_decomposition(0.10, method="fortran")

    # Newer method fails properly
    from capytaine.green_functions.abstract_green_function import GreenFunctionEvaluationError
    with pytest.raises((NotImplementedError, GreenFunctionEvaluationError)):
        gf.find_best_exponential_decomposition(0.10, method="python")
    with pytest.raises((NotImplementedError, GreenFunctionEvaluationError)):
        gf.find_best_exponential_decomposition(0.01, method="python")


@pytest.mark.xfail  # Find a way to pass parameters to find_best_exponential_decomposition to make sure they match
def test_python_and_fortran_prony_decomposition_for_green_function():
    gf = cpt.Delhommeau(finite_depth_prony_decomposition_method="python")
    decomp_default = gf.find_best_exponential_decomposition(1.0)
    decomp_f = gf.find_best_exponential_decomposition(1.0, method="fortran")
    decomp_p = gf.find_best_exponential_decomposition(1.0, method="python")
    assert np.allclose(decomp_default, decomp_p)
    assert np.allclose(decomp_p, decomp_f, rtol=0.2)


def test_prony_decomposition_specific_value():
    gf = cpt.Delhommeau(finite_depth_prony_decomposition_method="python")
    gf.find_best_exponential_decomposition(20.0)


def test_failure_unknown_method():
    gf = cpt.Delhommeau()
    with pytest.raises(ValueError):
        gf.find_best_exponential_decomposition(1.0, method="potato")


def test_fingreen3D():
    k = 5.0
    h = 2.0
    mesh = cpt.mesh_sphere(radius=0.5, center=(0.0, 0.0, -1.0), resolution=(4, 4))
    S, D = cpt.FinGreen3D().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    S_ref, D_ref = cpt.Delhommeau().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    np.testing.assert_allclose(S, S_ref, rtol=1e-1, atol=1e-1)
    np.testing.assert_allclose(D, D_ref, rtol=1e-1, atol=1e-1)


@pytest.mark.parametrize("omega", [0.1, 1.0, 10.0])
@pytest.mark.parametrize("depth", [0.1, 1.0, 10.0])
def test_fingreen3D_roots(omega, depth):
    from scipy.optimize import newton, brentq
    nk = 100
    omega2_over_g = omega**2/9.81
    all_roots = cpt.FinGreen3D().fortran_core.green_wave.compute_dispersion_roots(nk, omega2_over_g, depth)
    wavenumber = newton(lambda x: omega2_over_g - x*np.tanh(x*depth), omega2_over_g)
    assert np.isclose(wavenumber, all_roots[0], rtol=1e-3)

    roots_idempotent = newton(lambda x: -omega2_over_g - x*np.tan(x*depth), all_roots[1:])
    np.testing.assert_allclose(roots_idempotent, all_roots[1:], rtol=1e-3)

    def root(i_root):
        return brentq(lambda y: omega2_over_g*depth + y*np.tan(y), (2*i_root+1)*np.pi/2 + 1e-10, (2*i_root+2)*np.pi/2 - 1e-10)/depth
    roots_2 = [root(i_root) for i_root in range(nk-1)]
    np.testing.assert_allclose(roots_2, all_roots[1:], rtol=1e-4)

    all_roots_python = cpt.FinGreen3D().compute_dispersion_relation_roots(nk, wavenumber, depth)
    np.testing.assert_allclose(all_roots_python, all_roots, rtol=1e-3)
