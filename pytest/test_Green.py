#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

from itertools import product, combinations

import pytest

import numpy as np
from scipy.special import exp1
from scipy.optimize import newton

import capytaine.NemohCore as NemohCore


def E1(z):
    return np.exp(-z)*NemohCore.initialize_green_wave.gg(z)

@pytest.mark.parametrize("x", np.linspace(-20, -1, 4))
@pytest.mark.parametrize("y", np.linspace(-10, 10, 4))
def test_GG(x, y):
    # Test some properties of the function according to [Del, p.367].
    z = x + 1j*y
    assert np.isclose(E1(z), exp1(z), rtol=1e-3)

    # (A3.5)
    assert np.isclose(E1(x - 1j*y), np.conjugate(E1(x + 1j*y)), rtol=1e-3)

    # TODO: test if d/dz (e^z E1(z)) = e^z E1(z) - 1/z


@pytest.mark.parametrize("omega", np.linspace(0.1, 5.0, 2))
@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_green_function(omega, depth):
    g = 9.8

    if depth == np.infty:
        wavenumber = omega**2 / g
    else:
        wavenumber = newton(lambda x: x*np.tanh(x) - omega**2*depth/g, x0=1.0)/depth

    XR, XZ, APD = NemohCore.initialize_green_wave.initialize_tabulated_integrals(328, 46, 251)
    if depth < np.infty:
        ambda, ar, nexp = NemohCore.old_prony_decomposition.lisc(omega**2 * depth/g, wavenumber * depth)

    def g(w, Xi, Xj):
        if depth == np.infty:
            return NemohCore.green_wave.wave_part_infinite_depth(w, Xi, Xj, XR, XZ, APD)[0]
        else:
            return NemohCore.green_wave.wave_part_finite_depth(w, Xi, Xj, depth, XR, XZ, APD, ambda, ar, 31)[0]

    def dg(w, Xi, Xj):
        return NemohCore.green_wave.wave_part_infinite_depth(w, Xi, Xj, XR, XZ, APD)[1]

    def reflect(X):
        Y = X.copy()
        Y[0:2] = -Y[0:2]
        return Y

    x_range = np.linspace(-2.0, 2.0, 3)
    y_range = np.linspace(-2.0, 2.0, 3)
    z_range = np.linspace(-9.5, -0.5, 4)
    for Xi, Xj in combinations(product(x_range, y_range, z_range), 2):
        assert np.isclose(g(wavenumber, np.array(Xi), np.array(Xj)),
                          g(wavenumber, np.array(Xj), np.array(Xi)),
                          rtol=1e-4)
        assert np.allclose(dg(wavenumber, np.array(Xi), np.array(Xj)),
                           reflect(dg(wavenumber, np.array(Xj), np.array(Xi))),
                           rtol=1e-4)


