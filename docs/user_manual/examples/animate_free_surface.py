#!/usr/bin/env python

import logging

from capytaine import Nemoh, DiffractionProblem, RadiationProblem
from capytaine.geometric_bodies import Sphere, FreeSurface
from capytaine.tools.vtk.free_surface_animation import Animation

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

# Generate the mesh of a sphere
full_sphere = Sphere(
    radius=5, center=(0, 0, -2.0),  # Size and positions
    ntheta=30, nphi=30,             # Fineness of the mesh
)

# Keep only the immersed part of the mesh
sphere = full_sphere.get_immersed_part()

# Set up and solve problem
problem = DiffractionProblem(body=sphere, angle=0.0, omega=2.0)
solver = Nemoh()
result = solver.solve(problem, keep_details=True)

# Define a mesh of the free surface and compute the free surface elevation
fs_mesh = FreeSurface(x_range=(-50, 50), y_range=(-50, 50), nx=100, ny=100)
fs_elevation = solver.get_free_surface_elevation(result, fs_mesh)

# Add incoming waves
fs_elevation = fs_elevation + fs_mesh.incoming_waves(result)

# Run the animation
animation = Animation(
    result,
    fs_mesh,
    fs_mesh.elevation_at_nodes(fs_elevation),
    display_body=full_sphere,
)
animation.run()

