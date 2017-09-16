#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

import capytaine._Green as _G
from capytaine._Wavenumber import invert_xtanhx

mesh_resolution = 100

# depth = np.infty
depth = 10

omega = 1.0
g = 9.81

#################
#  Computation  #
#################

_G.initialize_green_2.initialize_green()
if depth == np.infty:
    wavenumber = omega**2/g
else:
    wavenumber = invert_xtanhx(omega**2*depth/g)/depth
    _G.initialize_green_2.lisc(omega**2*depth/g, wavenumber*depth)

source = np.asarray([0.0,  0.0, -5.0])

X_range = np.linspace(-2/wavenumber, 2/wavenumber, mesh_resolution)
if depth == np.infty:
    Z_range = np.linspace(-10.0, 0.1, mesh_resolution)
else:
    Z_range = np.linspace(-depth-0.2, 0.2, mesh_resolution)
X_mesh, Z_mesh = np.meshgrid(X_range, Z_range)

green = np.zeros((len(X_range), len(Z_range)), dtype=np.complex64)

for i, x in enumerate(X_range):
    for j, z in enumerate(Z_range):

        p = np.array([x, 0.0, z])

        green[i, j] += 1/(norm(source-p))

        if depth == np.infty:
            p_mirror = np.array([x, 0.0, -z])
            green[i, j] += 1/(norm(source-p_mirror))

            green[i, j] += _G.green_2.vnsinfd(wavenumber, source, p, 1.0)[0]

        else:
            p_mirror = np.array([x, 0.0, -2*depth-z])
            green[i, j] += 1/(norm(source-p_mirror))

            green[i, j] += _G.green_2.vnsfd(wavenumber, source, p, 1.0, depth)[0]

##########
#  Plot  #
##########

plt.figure()
plt.contourf(
    X_mesh,
    Z_mesh,
    np.log(np.real(green.T)),
    10
)
plt.colorbar()
plt.title("Green function")

plt.figure()
CS = plt.contour(
    X_mesh,
    Z_mesh,
    np.gradient(np.real(green.T), axis=0),
    levels=np.arange(-0.05, 0.05, 0.01),
)
plt.clabel(CS)
plt.title("Green function gradient w.r.t. $z$")

plt.figure()
CS = plt.contour(
    X_mesh,
    Z_mesh,
    omega**2/g * np.real(green.T) - np.gradient(np.real(green.T), axis=0),
    levels=np.arange(-0.05, 0.05, 0.01),
)
plt.clabel(CS)
plt.title("Green function, free surface boundary condition")

plt.show()
