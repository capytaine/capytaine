#!/usr/bin/env python
# coding: utf-8
"""
Example computation: plot the free surface for wave diffraction around a sphere
"""

import logging

from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.free_surface import FreeSurface

from capytaine.problems import DiffractionProblem
from capytaine.Nemoh import Nemoh
from capytaine.tools.Airy_wave import Airy_wave_potential
from capytaine.tools.VTK_free_surface_animation import Animation

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s:\t%(message)s", datefmt="%H:%M:%S")

# Initialize mesh and solver
sphere = Sphere(radius=5, ntheta=40, nphi=40, center=(0, 0, -2.0), clip_free_surface=True)
solver = Nemoh()

# Solve diffraction problem
problem = DiffractionProblem(body=sphere, angle=0.0, omega=2.0)
result = solver.solve(problem, keep_details=True)

# Compute free surface elevation
fs = FreeSurface(x_range=(-50, 50), y_range=(-50, 50), nx=100, ny=100)
fs_elevation = solver.get_free_surface_elevation(result, fs, keep_details=True)

# Add incoming waves
fs_elevation = fs_elevation + 1j * problem.omega / problem.g * Airy_wave_potential(fs.mesh.faces_centers, result)

Animation(result, fs, fs.elevation_at_nodes(fs_elevation)).run()
