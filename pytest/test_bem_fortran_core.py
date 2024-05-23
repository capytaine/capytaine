"""Tests for the computation of the Green function using Fortran routines."""

import pytest
import numpy as np
import capytaine as cpt
from pytest import approx

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

def test_infinite_depth_gf():
    from scipy.special import j0
    gf = cpt.Delhommeau()
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                args[3])
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
    gf = cpt.Delhommeau(gf_singularities=gf_singularities)
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                gf.gf_singularities_index)[0]
    def nabla_wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                gf.gf_singularities_index)[1]
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


gfs = [
        cpt.Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        cpt.XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        ]


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_infinite_depth(gf):
    k = 1.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xi, xj, k,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals, gf.gf_singularities_index
    )
    g2, dg2 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xj, xi, k,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals, gf.gf_singularities_index
    )
    assert g1 == approx(g2)
    assert dg1[0:2] == approx(-dg2[0:2])
    assert dg1[2] == approx(dg2[2])


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_finite_depth_no_prony(gf):
    k = 1.0
    depth = 5.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xi, xj, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        np.zeros(1), np.zeros(1), 1
    )
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xj, xi, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        np.zeros(1), np.zeros(1), 1
    )
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_finite_depth(gf):
    k = 1.0
    depth = 10.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    ambda, a, nexp = gf.fortran_core.old_prony_decomposition.lisc(k*depth*np.tanh(k*depth), k*depth)
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xi, xj, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        ambda, a, 31
    )
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xj, xi, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        ambda, a, 31
    )
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


def test_floating_point_precision():
    assert cpt.Delhommeau(floating_point_precision="float64").tabulated_integrals.dtype == np.float64
    assert cpt.Delhommeau(floating_point_precision="float32").tabulated_integrals.dtype == np.float32
