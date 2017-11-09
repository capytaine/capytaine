#!/usr/bin/env python
# coding: utf-8
"""
Example computation: plot the free surface for wave diffraction around a sphere
"""

import logging

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from capytaine.reference_bodies import generate_sphere, generate_free_surface
from capytaine.problems import DiffractionProblem
from capytaine.Nemoh import Nemoh
from matplotlib.patches import Circle

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s:\t%(message)s")

# Initialize mesh and solver
sphere = generate_sphere(radius=5, ntheta=40, nphi=40, z0=-1.0, clip_free_surface=True)
solver = Nemoh()

# Solve diffraction problem
problem = DiffractionProblem(body=sphere, angle=0.0, omega=2.0)
results = solver.solve(problem, keep_details=True)

# Compute free surface elevation
fs_mesh = generate_free_surface(width=100.0, length=100.0, nw=100, nl=100)
fs = solver.get_free_surface(problem, fs_mesh)

# Add incoming waves
fs = fs + 1j*problem.omega/problem.g * problem.Airy_wave_potential(fs_mesh.faces_centers)

# Plot free surface elevation
X = fs_mesh.faces_centers[:, 0].reshape(100, 100)
Y = fs_mesh.faces_centers[:, 1].reshape(100, 100)
fs = fs.reshape(100, 100)
scale = np.abs(fs).max()

fig = plt.figure()
ax = plt.gca()

nbi = 40 # Number of images in the animation

def animate(i):
    ax.clear()

    # Draw body shape
    body = Circle((0,0), 5, facecolor="w", edgecolor='w')
    ax.add_artist(body)

    # Draw free surface elevation
    ax.contourf(X, Y, np.abs(fs)*np.cos(np.angle(fs) - 2*np.pi*i/nbi), 30, vmin=-scale, vmax=scale)
    return ax

plt.contourf(X, Y, np.real(fs), 100, vmin=-scale, vmax=scale)
plt.colorbar()

plt.axis('equal')
ani = animation.FuncAnimation(fig, animate, nbi, interval=100/25, blit=False)

plt.show()

