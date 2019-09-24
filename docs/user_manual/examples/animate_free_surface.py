#!/usr/bin/env python

import logging

import capytaine as cpt
from capytaine.ui.vtk.animation import Animation

logging.basicConfig(level=logging.INFO, format="%(levelname)s:\t%(message)s")

# Generate the mesh of a sphere
full_sphere = cpt.Sphere(
    radius=3, center=(0, 0, 0),  # Size and positions
    ntheta=20, nphi=20,          # Fineness of the mesh
)
full_sphere.add_translation_dof(name="Heave")

# Keep only the immersed part of the mesh
sphere = full_sphere.keep_immersed_part(inplace=False)
sphere.add_translation_dof(name="Heave")

# Set up and solve problem
solver = cpt.BEMSolver()

diffraction_problem = cpt.DiffractionProblem(body=sphere, wave_direction=0.0, omega=2.0)
diffraction_result = solver.solve(diffraction_problem)

radiation_problem = cpt.RadiationProblem(body=sphere, radiating_dof="Heave", omega=2.0)
radiation_result = solver.solve(radiation_problem)

# Define a mesh of the free surface and compute the free surface elevation
free_surface = cpt.FreeSurface(x_range=(-50, 50), y_range=(-50, 50), nx=150, ny=150)
diffraction_elevation_at_faces = solver.get_free_surface_elevation(diffraction_result, free_surface)
radiation_elevation_at_faces = solver.get_free_surface_elevation(radiation_result, free_surface)

# Add incoming waves
diffraction_elevation_at_faces = diffraction_elevation_at_faces + free_surface.incoming_waves(diffraction_result)

# Run the animations
animation = Animation(loop_duration=diffraction_result.period)
animation.add_body(full_sphere, faces_motion=None)
animation.add_free_surface(free_surface, faces_elevation=0.5*diffraction_elevation_at_faces)
animation.run(camera_position=(-30, -30, 30))  # The camera is oriented towards (0, 0, 0) by default.
# animation.save("path/to/the/video/file.ogv", camera_position=(-30, -30, 30))

animation = Animation(loop_duration=radiation_result.period)
animation.add_body(full_sphere, faces_motion=full_sphere.dofs["Heave"])
animation.add_free_surface(free_surface, faces_elevation=3.0*radiation_elevation_at_faces)
animation.run(camera_position=(-30, -30, 30))
# animation.save("path/to/the/video/file.ogv", camera_position=(-30, -30, 30))

