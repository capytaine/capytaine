#!/usr/bin/env python
# coding: utf-8
"""
Example computation: plot the free surface for wave diffraction around a sphere
"""

import logging

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Circle

from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.free_surface import FreeSurface

from capytaine.problems import DiffractionProblem
from capytaine.Nemoh import Nemoh
from capytaine.tools.Airy_wave import Airy_wave_potential

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s:\t%(message)s", datefmt="%H:%M:%S")

# Initialize mesh and solver
sphere = Sphere(radius=5, ntheta=40, nphi=40, center=(0, 0, -1.0), clip_free_surface=True)
solver = Nemoh()

# Solve diffraction problem
problem = DiffractionProblem(body=sphere, angle=0.0, omega=2.0)
result = solver.solve(problem, keep_details=True)

# Compute free surface elevation
fs = FreeSurface(x_range=(-50, 50), y_range=(-50, 50), nx=100, ny=100)
fs_elevation = solver.get_free_surface_elevation(result, fs)

# Add incoming waves
fs_elevation = fs_elevation + 1j * problem.omega / problem.g * Airy_wave_potential(fs.mesh.faces_centers, result)

# Plot free surface elevation
X = fs.mesh.faces_centers[:, 0].reshape(fs.nx, fs.ny)
Y = fs.mesh.faces_centers[:, 1].reshape(fs.nx, fs.ny)
fs_elevation = fs_elevation.reshape(fs.nx, fs.ny)
scale = np.abs(fs_elevation).max()

fig = plt.figure()
ax = plt.gca()

nbi = 40 # Number of images in the animation


def animate(i):
    ax.clear()

    # Draw body shape
    body = Circle((0,0), 5, facecolor="w", edgecolor='w')
    ax.add_artist(body)

    # Draw free surface elevation
    ax.contourf(X, Y, np.abs(fs_elevation) * np.cos(np.angle(fs_elevation) - 2 * np.pi * i / nbi), 30, vmin=-scale, vmax=scale)
    return ax


plt.contourf(X, Y, np.real(fs_elevation), 100, vmin=-scale, vmax=scale)
plt.colorbar()

plt.axis('equal')
ani = animation.FuncAnimation(fig, animate, nbi, interval=100/25, blit=False)

plt.show()

