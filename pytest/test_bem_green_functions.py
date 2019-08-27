#!/usr/bin/env python
# coding: utf-8
"""Tests for the computation of the Green function using Fortran routines from Nemoh."""

from itertools import product, combinations

import pytest
from hypothesis import given, assume
from hypothesis.strategies import one_of, just, floats
from hypothesis.extra.numpy import arrays

import numpy as np
from scipy.special import exp1
from scipy.misc import derivative
from scipy.optimize import newton

from capytaine.green_functions import Delhommeau_f90, XieDelhommeau_f90


def E1(z):
    return np.exp(-z)*Delhommeau_f90.initialize_green_wave.gg(z)

@given(floats(min_value=-1e2, max_value=-1e-2), floats(min_value=-1e2, max_value=1e2))
def test_GG(x, y):
    z = x + 1j*y

    # Compare with Scipy implementation
    assert np.isclose(E1(z), exp1(z), rtol=1e-3)

    # Test property (A3.5) of the function according to [Del, p.367]. 
    if y != 0.0:
        assert np.isclose(E1(np.conjugate(z)), np.conjugate(E1(z)), rtol=1e-3)

    # BROKEN...
    # def derivative_of_complex_function(f, z, **kwargs):
    #     direction = 1
    #     return derivative(lambda eps: f(z + eps*direction), 0.0, **kwargs)
    # if abs(x) > 1 and abs(y) > 1:
    #     assert np.isclose(
    #         derivative_of_complex_function(Delhommeau_f90.initialize_green_wave.gg, z),
    #         Delhommeau_f90.initialize_green_wave.gg(z) - 1.0/z,
    #         atol=1e-2)


tabulation = {
    Delhommeau_f90: Delhommeau_f90.initialize_green_wave.initialize_tabulated_integrals(328, 46, 251),
    XieDelhommeau_f90: XieDelhommeau_f90.initialize_green_wave.initialize_tabulated_integrals(328, 46, 251),
}

def test_tabulations():
    assert np.allclose(tabulation[Delhommeau_f90][2][:, :, 0, 0], tabulation[XieDelhommeau_f90][2][:, :, 0, 0])
    assert np.allclose(tabulation[Delhommeau_f90][2][:, :, 1, 0], tabulation[XieDelhommeau_f90][2][:, :, 1, 0])
    assert not np.allclose(tabulation[Delhommeau_f90][2][:, :, 0, 1], tabulation[XieDelhommeau_f90][2][:, :, 0, 1])
    assert np.allclose(tabulation[Delhommeau_f90][2][:, :, 1, 1], tabulation[XieDelhommeau_f90][2][:, :, 1, 1])


points = arrays(np.float, (3,),
                elements=floats(min_value=-1e5, max_value=1e5, allow_infinity=False, allow_nan=False)
                ).filter(lambda x: x[2] < -1e-4)
cores = one_of(just(Delhommeau_f90), just(XieDelhommeau_f90))
frequencies = floats(min_value=1e-1, max_value=1e1)
depths = one_of(floats(min_value=10.0, max_value=100.0), just(np.infty))

gravity = 9.8

def wave_part_Green_function(Xi, Xj, omega, depth, core=Delhommeau_f90):
    if depth == np.infty:
        wavenumber = omega**2 / gravity
    else:
        wavenumber = newton(lambda x: x*np.tanh(x) - omega**2*depth/gravity, x0=1.0)/depth

    if depth < np.infty:
        ambda, ar, nexp = core.old_prony_decomposition.lisc(omega**2 * depth/gravity, wavenumber * depth)

    if depth == np.infty:
        return core.green_wave.wave_part_infinite_depth(wavenumber, Xi, Xj, *tabulation[core])
    else:
        return core.green_wave.wave_part_finite_depth(wavenumber, Xi, Xj, depth, *tabulation[core], ambda, ar, 31)

@given(points, points, frequencies, depths, cores)
def test_symmetry_of_the_Green_function(X1, X2, omega, depth, core):
    assert np.isclose(wave_part_Green_function(X1, X2, omega, depth, core=core)[0],
                      wave_part_Green_function(X2, X1, omega, depth, core=core)[0],
                      rtol=1e-4)

@given(points, points, frequencies, cores)
def test_symmetry_of_the_derivative_of_the_Green_function(X1, X2, omega, core):
    assert np.allclose( wave_part_Green_function(X1, X2, omega, np.infty, core=core)[1][0:2],
                       -wave_part_Green_function(X2, X1, omega, np.infty, core=core)[1][0:2],
                      rtol=1e-4)
    assert np.isclose(wave_part_Green_function(X1, X2, omega, np.infty, core=core)[1][2],
                      wave_part_Green_function(X2, X1, omega, np.infty, core=core)[1][2],
                      rtol=1e-4)


