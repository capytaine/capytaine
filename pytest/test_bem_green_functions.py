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

from capytaine import Delhommeau, XieDelhommeau


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
        Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251),
        XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251),
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


@given(points, points, frequencies)
def test_compare_output_of_Delhommeau_and_XieDelhommeau(X1, X2, omega):
    # Unchanged parts of the Green function:
    assert np.allclose(
            wave_part_Green_function(X1, X2, omega, np.infty, gfs[0])[1][0:2],
            wave_part_Green_function(X1, X2, omega, np.infty, gfs[1])[1][0:2],
            rtol=1e-5)

