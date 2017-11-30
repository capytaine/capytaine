#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

from itertools import product

import numpy as np

import capytaine._Green as _Green

def gg(z):
    return _Green.initialize_green_2.gg(z, 0j)

def test_GG():
    # Test some properties of the function according to [Del, p.367].

    # (A3.5)
    for x, y in product(np.linspace(-10, 10, 10), np.linspace(-10, 10, 10)):
        assert gg(x - 1j*y) == np.conjugate(gg(x + 1j*y))

    # TODO: test if d/dz (e^z gg(z)) = e^z gg(z) - 1/z
