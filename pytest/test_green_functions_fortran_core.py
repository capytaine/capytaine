"""Tests for the computation of the Green function using Fortran routines."""

import pytest

import numpy as np
import capytaine as cpt

# def test_exponential_integral(x, y):
#     # Compare with Scipy implementation
#     from scipy.special import exp1
#     z = x + 1j*y
#
#     gf = cpt.Delhommeau()
#     def E1(z):
#         return np.exp(-z)*gf.fortran_core.delhommeau_integrals.exp_e1(z)
#
#     assert np.isclose(E1(z), exp1(z), rtol=1e-3)

#     # Test property (A3.5) of the function according to [Del, p.367].
#     if y != 0.0:
#         assert np.isclose(E1(np.conjugate(z)), np.conjugate(E1(z)), rtol=1e-3)


def test_infinite_depth_gf_singularities_variants():
    """Test how the values of the variants of the infinite depth Green function are related."""
    from scipy.special import j0
    gf = cpt.Delhommeau()
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], *gf.all_tabulation_parameters, args[3])
    xi = np.array([0.0, 0.0, -0.5])
    xj = np.array([0.0, 1.0, -0.5])
    k = 2.0
    g_lf, nablag_lf = wave_part(xi, xj, k, 0)
    g_hf, nablag_hf = wave_part(xi, xj, k, 1)
    r = k * np.linalg.norm(xi[:2] - xj[:2])
    z = k * (xi[2] + xj[2])
    assert g_hf.real + 2*k/np.hypot(r, z) == pytest.approx(g_lf.real, rel=1e-4)
    assert g_lf.imag == pytest.approx(g_hf.imag)
    assert g_lf.imag == pytest.approx(2*np.pi*k*np.exp(z)*j0(r), rel=1e-4)


@pytest.mark.parametrize("gf_singularities", ["low_freq", "high_freq"])
def test_infinite_depth_gf_finite_differences(gf_singularities):
    """Compare the gradient of the Green function as returned by the code with the finite difference gradient"""
    # TODO: same in finite depth
    gf = cpt.Delhommeau(gf_singularities=gf_singularities)
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], *gf.all_tabulation_parameters, gf.gf_singularities_fortran_enum[gf.gf_singularities])[0]
    def nabla_wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], *gf.all_tabulation_parameters, gf.gf_singularities_fortran_enum[gf.gf_singularities])[1]
    k = 1.0
    xi = np.array([0.0, 0.0, -0.5])
    xj = np.array([1.0, 1.0, -0.5])
    x_eps = np.array([1e-3, 0.0, 0.0])
    y_eps = np.array([0.0, 1e-3, 0.0])
    z_eps = np.array([0.0, 0.0, 1e-3])
    g = wave_part(xi, xj, k)
    g_x = (wave_part(xi + x_eps, xj, k) - g)/np.linalg.norm(x_eps)
    g_y = (wave_part(xi + y_eps, xj, k) - g)/np.linalg.norm(y_eps)
    g_z = (wave_part(xi + z_eps, xj, k) - g)/np.linalg.norm(z_eps)
    nabla_g = np.array([g_x, g_y, g_z])
    nablag_ref = nabla_wave_part(xi, xj, k)
    assert nablag_ref == pytest.approx(nabla_g, abs=1e-2)


def test_value_of_inf_depth_green_function_wave_part():
    gf = cpt.Delhommeau()
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], *gf.all_tabulation_parameters, gf.fortran_core.constants.low_freq)
    k = 1.0
    xi = np.array([0.0, 0.0, -0.5])
    xj = np.array([1.0, 1.0, -0.5])
    g, nabla_g = wave_part(xi, xj, k)
    assert g == pytest.approx(-2.04 + 1.29j, abs=0.1)
    assert nabla_g == pytest.approx([0.25+0.88j,  0.25+0.88j, -0.88+1.29j], abs=0.1)


gfs = [
        cpt.Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        cpt.XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        # TODO: add more relevant Green functions
        ]

@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_infinite_depth(gf):
    k = 1.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xi, xj, k, *gf.all_tabulation_parameters, gf.gf_singularities_fortran_enum[gf.gf_singularities]
    )
    g2, dg2 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xj, xi, k, *gf.all_tabulation_parameters, gf.gf_singularities_fortran_enum[gf.gf_singularities]
    )
    assert g1 == pytest.approx(g2)
    assert dg1[0:2] == pytest.approx(-dg2[0:2])
    assert dg1[2] == pytest.approx(dg2[2])


def test_floating_point_precision():
    assert cpt.Delhommeau(floating_point_precision="float64").tabulated_integrals.dtype == np.float64
    assert cpt.Delhommeau(floating_point_precision="float32").tabulated_integrals.dtype == np.float32


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
    rankine_g_once = gf.evaluate_rankine_only(center.reshape(1, -1), mesh)[0]
    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(11, 11))
    # Use odd number of panels to avoid having the center on a corner of a panel, which is not defined for the strong singularity of the derivative.
    rankine_g_parts = gf.evaluate_rankine_only(center.reshape(1, -1), mesh)[0]
    assert rankine_g_once == pytest.approx(rankine_g_parts.sum(), abs=1e-2)


def test_rankine_term_alone_or_in_zero_freq():
    import capytaine as cpt
    gf = cpt.Delhommeau()
    collocation_point = np.array([1.0, 0.0, -1.0]).reshape(1, -1)
    mirror_collocation_point = np.array([1.0, 0.0, 1.0]).reshape(1, -1)
    mesh = cpt.mesh_rectangle(size=(1.0, 1.0), center=(0.0, 0.0, -1.0), resolution=(1, 1))
    rankine_only = gf.evaluate_rankine_only(collocation_point, mesh)[0] + gf.evaluate_rankine_only(mirror_collocation_point, mesh)[0]
    zero_freq = gf.evaluate(collocation_point, mesh, water_depth=np.inf, wavenumber=0.0)[0]
    assert rankine_only == pytest.approx(zero_freq)


def test_rankine_term_alone_or_in_inf_freq():
    import capytaine as cpt
    gf = cpt.Delhommeau()
    collocation_point = np.array([1.0, 0.0, -1.0]).reshape(1, -1)
    mirror_collocation_point = np.array([1.0, 0.0, 1.0]).reshape(1, -1)
    mesh = cpt.mesh_rectangle(size=(1.0, 1.0), center=(0.0, 0.0, -1.0), resolution=(1, 1))
    rankine_only = gf.evaluate_rankine_only(collocation_point, mesh)[0] - gf.evaluate_rankine_only(mirror_collocation_point, mesh)[0]
    inf_freq = gf.evaluate(collocation_point, mesh, water_depth=np.inf, wavenumber=np.inf)[0]
    assert rankine_only == pytest.approx(inf_freq)


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
    # TODO: also for not adjoint double layer
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


def test_liangwunoblesse_wave_term():
    core = cpt.LiangWuNoblesseGF().fortran_core
    gf, gf_r = core.liangwunoblessewaveterm.havelockgf(1.0, -1.0)
    re_gf, _, im_gf, re_gf_r, im_gf_r = core.delhommeau_integrals(1.0, -1.0, 1000)
    assert gf == pytest.approx(- re_gf - im_gf*1j, rel=5e-2)
    assert gf_r == pytest.approx(- re_gf_r - im_gf_r*1j, rel=5e-2)


@pytest.mark.parametrize("wavenumber", [0.0, 1.0, np.inf])
def test_full_liangwunoblesse(wavenumber):
    gf = cpt.LiangWuNoblesseGF()
    ref_gf = cpt.Delhommeau()
    mesh = cpt.mesh_sphere().immersed_part()
    S, K = gf.evaluate(mesh, mesh, 0.0, np.inf, wavenumber)
    ref_S, ref_K = ref_gf.evaluate(mesh, mesh, 0.0, np.inf, wavenumber)
    np.testing.assert_allclose(S, ref_S, rtol=1e-2)
    np.testing.assert_allclose(K, ref_K, rtol=1e-2)
