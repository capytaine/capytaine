#!/usr/bin/env python
# coding: utf-8
"""
Example computation: plot the free surface elevation
"""

import sys
import logging

from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.free_surface import FreeSurface

from capytaine.problems import DiffractionProblem, RadiationProblem
from capytaine.Nemoh import Nemoh

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s:\t%(message)s", datefmt="%H:%M:%S")

# Initialize mesh and solver
sphere = Sphere(radius=5, ntheta=40, nphi=40, center=(0, 0, -2.0), clip_free_surface=True)
sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")

problem = RadiationProblem(body=sphere, radiating_dof="Surge", omega=2.0)
# problem = DiffractionProblem(body=sphere, angle=0.0, omega=2.0)

# Solve problem
solver = Nemoh()
result = solver.solve(problem, keep_details=True)

# Compute free surface elevation
fs = FreeSurface(x_range=(-50, 50), y_range=(-50, 50), nx=40, ny=40)
fs_elevation = solver.get_free_surface_elevation(result, fs)

# Add incoming waves
if isinstance(problem, DiffractionProblem):
    fs_elevation = fs_elevation + fs.incoming_waves(result)

# Animation
if "vtk" in sys.modules:
    from capytaine.tools.VTK_free_surface_animation import Animation
    Animation(result, fs, fs.elevation_at_nodes(fs_elevation)).run()

else:
    from matplotlib.patches import Circle
    from capytaine.tools.mpl_free_surface_animation import animation_matplotlib
    animation_matplotlib(result, fs, fs_elevation,
                         body_shape=Circle((0,0), 5, facecolor="w", edgecolor='w'))

