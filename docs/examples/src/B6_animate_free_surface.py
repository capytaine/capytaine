"""This example shows how to animate the free surface elevation and pressure on a floating body using vedo."""

import numpy as np
import capytaine as cpt
from capytaine.bem.airy_waves import (
    airy_waves_free_surface_elevation,
    airy_waves_pressure,
)
from capytaine.ui.vedo_animations import Animation


# Parameters
omega = 3.0
T = 2 * np.pi / omega


## RESOLUTION WITH CAPYTAINE
positions = [
    (np.cos(angle) * 5.0, np.sin(angle) * 5.0, 0.0)
    for angle in [0.0, 2 * np.pi / 3, -2 * np.pi / 3]
]
meshes = [
    cpt.mesh_vertical_cylinder(length=4.0, radius=2.0, center=p, faces_max_radius=0.3)
    for p in positions
]
mesh = cpt.Mesh.join_meshes(*meshes)
sphere = cpt.FloatingBody(mesh=mesh.immersed_part(), lid_mesh=mesh.generate_lid())
solver = cpt.BEMSolver()
diffraction_problem = cpt.DiffractionProblem(body=sphere, wave_direction=0.0, omega=omega)
print("Solving the diffraction problem...")
diffraction_result = solver.solve(diffraction_problem)
print("Diffraction problem solved.")


## POST-PROCESSING THE PRESSURE AND FREE SURFACE ELEVATION
free_surface_mesh = cpt.mesh_rectangle(
    size=(100.0, 100.0), center=(0.0, 0.0, 0.0), resolution=(300, 300)
)
amplitude = 0.2
print("Computing the free surface elevation and pressure...")
diffraction_elevation_at_faces = amplitude * (
    solver.compute_free_surface_elevation(free_surface_mesh.vertices, diffraction_result)
    + airy_waves_free_surface_elevation(free_surface_mesh.vertices, diffraction_problem)
)
pressure_on_sphere = amplitude * (
    solver.compute_pressure(sphere.mesh.vertices, diffraction_result)
    + airy_waves_pressure(sphere.mesh.vertices, diffraction_problem)
)
print("Free surface elevation and pressure computed.")


## ANIMATION WITH VEDO
anim = Animation(loop_duration=T, fps=24)
anim.add_body(sphere, vertices_pressure=pressure_on_sphere)
anim.add_free_surface(free_surface_mesh, vertices_elevation=diffraction_elevation_at_faces)
anim.save("diffraction.mp4", size=(1024, 512))
