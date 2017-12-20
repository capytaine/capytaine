#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

from itertools import product, combinations

import pytest
import numpy as np

import capytaine._Green as _Green
from capytaine._Wavenumber import invert_xtanhx


def test_GG():
    # Test some properties of the function according to [Del, p.367].

    def E1(z):
        return _Green.initialize_green_2.gg(z, 0j)

    # (A3.5)
    for x, y in product(np.linspace(-10, 10, 10), np.linspace(-10, 10, 10)):
        assert E1(x - 1j*y) == np.conjugate(E1(x + 1j*y))

    # TODO: test if d/dz (e^z E1(z)) = e^z E1(z) - 1/z


@pytest.mark.parametrize("omega", np.linspace(0.1, 5.0, 2))
@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_green_function(omega, depth):
    g = 9.8

    if depth == np.infty:
        wavenumber = omega**2 / g
    else:
        wavenumber = invert_xtanhx(omega**2 * depth/g) / depth

    XR = _Green.initialize_green_2.initialize_green()
    if depth < np.infty:
        _Green.initialize_green_2.lisc(omega**2 * depth/g, wavenumber * depth)

    def g(w, Xi, Xj):
        if depth == np.infty:
            return _Green.green_2.vnsinfd(w, Xi, Xj, XR)[0]
        else:
            return _Green.green_2.vnsfd(w, Xi, Xj, depth, XR)[0]

    def dg(w, Xi, Xj):
        return _Green.green_2.vnsinfd(w, Xi, Xj, XR)[1]

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


