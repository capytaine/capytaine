#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

import capytaine._Green as _G
from capytaine._Wavenumber import invert_xtanhx

mesh_resolution = 200

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

X_range = np.linspace(-4/wavenumber, 4/wavenumber, 2*mesh_resolution)
if depth == np.infty:
    Z_range = np.linspace(-10.0, 0.1, mesh_resolution)
else:
    Z_range = np.linspace(-depth-0.2, 0.2, mesh_resolution)
X_mesh, Z_mesh = np.meshgrid(X_range, Z_range)

green = np.zeros((len(Z_range), len(X_range)), dtype=np.complex64)

for i, x in enumerate(X_range):
    for j, z in enumerate(Z_range):

        p = np.array([x, 0.0, z])
        # green[j, i] += 1/(4*np.pi*norm(source-p))

        if depth == np.infty:
            p_mirror = np.array([x, 0.0, -z])
            # green[j, i] += 1/(4*np.pi*norm(source-p_mirror))
            green[j, i] += _G.green_2.vnsinfd(wavenumber, source, p)[0]

        else:
            p_mirror = np.array([x, 0.0, -2*depth-z])
            # green[j, i] += 1/(4*np.pi*norm(source-p_mirror))
            green[j, i] += _G.green_2.vnsfd(wavenumber, source, p, depth)[0]

##########
#  Plot  #
##########

plt.figure()
plt.contourf(X_mesh, Z_mesh, np.abs(green), 30)
plt.colorbar()
plt.title("Green function")

# plt.figure()
# laplace = np.gradient(np.gradient(np.real(green), axis=0), axis=0) \
#     + 2*np.gradient(np.gradient(np.real(green), axis=1), axis=1)
# plt.contourf(X_mesh, Z_mesh, laplace, 10)
# plt.colorbar()
# plt.title("Laplacian of Green function")

if depth < np.infty:
    plt.figure()
    indices = np.where(Z_range < -depth+0.2)[0]
    for i, index in enumerate(indices):
        plt.plot(
            X_range,
            np.gradient(np.real(green), axis=0)[index, :],
            color=(i/len(indices), 0, 0),
            label=f"z={Z_range[index]:.2f}",
        )
    plt.legend()
    plt.title("Gradient w.r.t. $z$")

plt.figure()
indices = np.where(Z_range > -0.2)[0]
BC = omega**2/g * np.real(green) - np.gradient(np.real(green), axis=0)
for i, index in enumerate(indices):
    plt.plot(X_range, BC[index, :], color=(0, i/len(indices), 0), label=f"z={Z_range[index]:.2f}")
plt.title("Boundary condition at the free surface")
plt.legend()

plt.show()
