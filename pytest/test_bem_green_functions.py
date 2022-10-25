#!/usr/bin/env python
# coding: utf-8
"""Tests for the computation of the Green function using Fortran routines from Nemoh."""

from itertools import product, combinations

import pytest
from hypothesis import given, assume
from hypothesis.strategies import one_of, just, floats
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import pi
from numpy.linalg import norm
from scipy.integrate import quad
from scipy.misc import derivative
from scipy.optimize import newton
from scipy.special import exp1

import capytaine as cpt


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
        cpt.Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251),
        cpt.XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251),
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
    Del = gfs[0].tabulated_integrals[:, :, 0, 1] - pi/R1
    Xie = gfs[1].tabulated_integrals[:, :, 0, 1]
    assert np.allclose(Del[abs(z) > 1], Xie[abs(z) > 1], atol=1e-3)

points = arrays(float, (3,), elements=floats(min_value=-10.0, max_value=-1e-2, allow_infinity=False, allow_nan=False))
methods = one_of(just(gfs[0]), just(gfs[1]))
frequencies = floats(min_value=1e-2, max_value=1e2)
depths = one_of(floats(min_value=10.0, max_value=100.0), just(np.infty))

gravity = 9.8

def wave_part_Green_function(Xi, Xj, omega, depth, method):
    if depth == np.infty:
        wavenumber = omega**2 / gravity
        return method.fortran_core.green_wave.wave_part_infinite_depth(Xi, Xj, wavenumber, method.tabulated_r_range, method.tabulated_z_range, method.tabulated_integrals)
    else:
        wavenumber = newton(lambda x: x*np.tanh(x) - omega**2*depth/gravity, x0=1.0)/depth
        ambda, ar, nexp = method.fortran_core.old_prony_decomposition.lisc(omega**2 * depth/gravity, wavenumber * depth)
        return method.fortran_core.green_wave.wave_part_finite_depth(Xi, Xj, wavenumber, depth, method.tabulated_r_range, method.tabulated_z_range, method.tabulated_integrals, ambda, ar, 31)


@given(points, points, frequencies, depths, methods)
def test_symmetry_of_the_Green_function(X1, X2, omega, depth, method):
    assert np.isclose(wave_part_Green_function(X1, X2, omega, depth, method)[0],
                      wave_part_Green_function(X2, X1, omega, depth, method)[0],
                      rtol=1e-4)

@given(points, points, frequencies, methods)
def test_symmetry_of_the_derivative_of_the_Green_function(X1, X2, omega, method):
    assert np.allclose( wave_part_Green_function(X1, X2, omega, np.infty, method)[1][0:2],
                       -wave_part_Green_function(X2, X1, omega, np.infty, method)[1][0:2],
                      rtol=1e-4)
    assert np.isclose(wave_part_Green_function(X1, X2, omega, np.infty, method)[1][2],
                      wave_part_Green_function(X2, X1, omega, np.infty, method)[1][2],
                      rtol=1e-4)


def test_floating_point_precision():
    assert cpt.Delhommeau(floating_point_precision="float64").tabulated_integrals.dtype == np.float64
    assert cpt.Delhommeau(floating_point_precision="float32").tabulated_integrals.dtype == np.float32


def test_no_tabulation():
    mesh = cpt.Sphere().mesh.keep_immersed_part()
    tabed_gf = cpt.Delhommeau()
    untabed_gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0)
    assert np.allclose(untabed_gf.evaluate(mesh, mesh)[0], tabed_gf.evaluate(mesh, mesh)[0], atol=1e-2)



def test_liang_wu_noblesse():
    import capytaine as cpt
    from capytaine.green_functions.delhommeau import LiangWuNoblesse
    delh = cpt.XieDelhommeau(tabulation_nr=100, tabulation_nz=100)
    tab = (delh.tabulated_r_range, delh.tabulated_z_range, delh.tabulated_integrals)
    lwn = LiangWuNoblesse()

    from itertools import product
    wavenumber = 1.0
    r_range = np.linspace(0.1, 100.0, 5); z_range = np.linspace(-100.0, -0.1, 5)
    for (r, z) in product(r_range, z_range):
        x, xi = np.array([r/np.sqrt(2), r/np.sqrt(2), z/2]), np.array([0.0, 0.0, z/2])
        dgf = delh.fortran_core.green_wave.collect_delhommeau_integrals(x, xi, wavenumber, *tab)
        lwngf = lwn.fortran_core.green_wave.collect_delhommeau_integrals(x, xi, wavenumber, *tab)
        print(f"{r:8.2f}, {z:8.2f}, {np.abs(dgf[0]-lwngf[0])/np.abs(dgf[0]):.4f}, {np.abs(dgf[1]-lwngf[1])/np.abs(dgf[1])}")

    m = cpt.Sphere(center=(0.0, 0.0, -2.0)).mesh
    S, K = LiangWuNoblesse().evaluate(m, m, 0.0, -np.infty, 1.0)
    S1, K1 = cpt.XieDelhommeau().evaluate(m, m, 0.0, -np.infty, 1.0)
    print(np.max(np.abs(S-S1)))
    print(np.max(np.abs(K-K1)))

test_liang_wu_noblesse()
