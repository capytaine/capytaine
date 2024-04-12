#!/usr/bin/env python
# coding: utf-8
"""Tests for the computation of the Green function using Fortran routines."""

import pytest
from pytest import approx

import numpy as np
import capytaine as cpt



# from numpy import pi
# from scipy.integrate import quad
# from scipy.misc import derivative
# from scipy.optimize import newton
# from scipy.special import exp1
# from hypothesis import given

# def E1(z):
#     return np.exp(-z)*Delhommeau_f90.initialize_green_wave.exp_e1(z)
#
# @given(floats(min_value=-1e2, max_value=-1e-2), floats(min_value=-1e2, max_value=1e2))
# def test_exponential_integral(x, y):
#     z = x + 1j*y
#
#     # Compare with Scipy implementation
#     assert np.isclose(E1(z), exp1(z), rtol=1e-3)
#
#     # Test property (A3.5) of the function according to [Del, p.367].
#     if y != 0.0:
#         assert np.isclose(E1(np.conjugate(z)), np.conjugate(E1(z)), rtol=1e-3)
#
#     # BROKEN...
#     # def derivative_of_complex_function(f, z, **kwargs):
#     #     direction = 1
#     #     return derivative(lambda eps: f(z + eps*direction), 0.0, **kwargs)
#     # if abs(x) > 1 and abs(y) > 1:
#     #     assert np.isclose(
#     #         derivative_of_complex_function(Delhommeau_f90.initialize_green_wave.gg, z),
#     #         Delhommeau_f90.initialize_green_wave.gg(z) - 1.0/z,
#     #         atol=1e-2)
#

# CHECK OF THE THEORY, NOT THE CODE ITSELF
# Not necessary to run every time...
#
# @given(r=floats(min_value=0, max_value=1e2),
#        z=floats(min_value=-16, max_value=-1),
#        k=floats(min_value=0.1, max_value=10)
#        )
# def test_zeta(r, z, k):
#     """Check expression k/π ∫ 1/ζ(θ) dθ = - 1/r """

#     def one_over_zeta(theta):
#         return np.real(1/(k*(z + 1j*r*np.cos(theta))))

#     int_one_over_zeta, *_ = quad(one_over_zeta, -pi/2, pi/2)

#     assert np.isclose(k/pi * int_one_over_zeta,
#                       -1/(np.sqrt(r**2 + z**2)),
#                       rtol=1e-4)


gfs = [
        cpt.Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        cpt.XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        ]

def test_compare_tabulations_of_Delhommeau_and_XieDelhommeau():
    assert np.allclose(gfs[0].tabulated_integrals[:, :, 0, 0],
                       gfs[1].tabulated_integrals[:, :, 0, 0])
    assert np.allclose(gfs[0].tabulated_integrals[:, :, 1, 0],
                       gfs[1].tabulated_integrals[:, :, 1, 0])
    assert np.allclose(gfs[0].tabulated_integrals[:, :, 1, 1],
                       gfs[1].tabulated_integrals[:, :, 1, 1])

    r_range = gfs[0].tabulated_r_range
    z_range = gfs[0].tabulated_z_range

    # Make 2D arrays
    r = r_range[:, None] * np.ones_like(z_range)[None, :]
    z = np.ones_like(r_range)[:, None] * z_range[None, :]
    R1 = np.hypot(r, z)

    # Compare Z1, for great values of Z, where both methods are as accurate
    Del = gfs[0].tabulated_integrals[:, :, 0, 1] - 1/R1
    Xie = gfs[1].tabulated_integrals[:, :, 0, 1]
    assert np.allclose(Del[abs(z) > 1], Xie[abs(z) > 1], atol=1e-3)


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_infinite_depth(gf):
    k = 1.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_infinite_depth(xi, xj, k, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals)
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_infinite_depth(xj, xi, k, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals)
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_finite_depth_no_prony(gf):
    k = 1.0
    depth = 5.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(xi, xj, k, depth, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals, np.zeros(1), np.zeros(1), 1)
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(xj, xi, k, depth, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals, np.zeros(1), np.zeros(1), 1)
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
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(xi, xj, k, depth, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals, ambda, a, 31)
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(xj, xi, k, depth, gf.tabulation_grid_shape_index, gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals, ambda, a, 31)
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


def test_floating_point_precision():
    assert cpt.Delhommeau(floating_point_precision="float64").tabulated_integrals.dtype == np.float64
    assert cpt.Delhommeau(floating_point_precision="float32").tabulated_integrals.dtype == np.float32


def test_no_tabulation():
    mesh = cpt.Sphere().mesh.keep_immersed_part()
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


@pytest.mark.parametrize("gf", gfs)
def test_exact_integration_of_rankine_terms(gf):
    center = np.array([0.0, 0.0, 0.0])
    area = 1.0

    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(1, 1))
    rankine_g_once = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(11, 11))
    # Use odd number of panels to avoid having the center on a corner of a panel, which is not defined for the strong singularity of the derivative.
    rankine_g_parts = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    assert rankine_g_once == pytest.approx(rankine_g_parts.sum(), abs=1e-2)
