#!/usr/bin/env python
# coding: utf-8

import numpy as np

import capytaine.NemohCore._Green as _G
from capytaine.NemohCore._Wavenumber import invert_xtanhx

omega = 1.0
depth = 10.0
g = 9.81
wavenumber = invert_xtanhx(omega**2*depth/g)/depth

_G.initialize_green_2.initialize_green()
_G.initialize_green_2.lisc(omega**2*depth/g, wavenumber*depth)

x_i = np.asarray([1.0,  2.0, -3.0])
x_j = np.asarray([-5.0, 1.0, -1.0])

print(_G.green_2.vnsinfd(wavenumber, x_i, x_j, 1.0))
print(_G.green_2.vnsinfd(wavenumber, x_j, x_i, 1.0))

# _G.green_2.vnsfd(wavenumber, x_i, x_j, 1.0, depth)
# _G.green_2.vnsfd(wavenumber, x_j, x_i, 1.0, depth)
