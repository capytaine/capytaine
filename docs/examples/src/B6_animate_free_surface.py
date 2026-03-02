"""This example shows how to animate the free surface elevation and pressure on a floating body using vedo."""

import numpy as np
from vedo import Plotter, Mesh, Video

import capytaine as cpt
from capytaine.bem.airy_waves import (
    airy_waves_free_surface_elevation,
    airy_waves_pressure,
)


# Parameters
omega = 3.0
T = 2 * np.pi / omega


## RESOLUTION WITH CAPYTAINE
positions = [
    (np.cos(angle) * 5.0, np.sin(angle) * 5.0, 0.0)
    for angle in [0.0, 2 * np.pi / 3, -2 * np.pi / 3]
]
meshes = [
    cpt.mesh_vertical_cylinder(length=3.0, radius=2.0, center=p, faces_max_radius=0.3)
    for p in positions
]
mesh = cpt.Mesh.join_meshes(*meshes)
sphere = cpt.FloatingBody(
    mesh=mesh.immersed_part(),
    lid_mesh=mesh.generate_lid(),
)
solver = cpt.BEMSolver()
diffraction_problem = cpt.DiffractionProblem(
    body=sphere, wave_direction=0.0, omega=omega
)
print("Solving the diffraction problem...")
diffraction_result = solver.solve(diffraction_problem)
print("Diffraction problem solved.")


## POST-PROCESSING THE PRESSURE AND FREE SURFACE ELEVATION
free_surface_mesh = cpt.mesh_rectangle(
    size=(100.0, 100.0), center=(0.0, 0.0, 0.0), resolution=(300, 300)
)
amplitude = 0.2
print("Computing the free surface elevation and pressure on the sphere...")
diffraction_elevation_at_faces = amplitude * (
    solver.compute_free_surface_elevation(
        free_surface_mesh.vertices, diffraction_result
    )
    + airy_waves_free_surface_elevation(free_surface_mesh.vertices, diffraction_problem)
)
pressure_on_sphere = amplitude * (
    solver.compute_pressure(sphere.mesh.vertices, diffraction_result)
    + airy_waves_pressure(sphere.mesh.vertices, diffraction_problem)
)
print("Free surface elevation and pressure computed.")


## ANIMATION WITH VEDO
cam = dict(
    pos=(-20.0, +20.0, -20.0),
    viewup=(0, 0, 1),
    focal_point=(0, 0, 0),
)
plt = Plotter(axes=1, size=(1024, 512), interactive=False)
vedo_sphere = Mesh([sphere.mesh.vertices, sphere.mesh.faces])
vedo_fs = Mesh([free_surface_mesh.vertices, free_surface_mesh.faces])
plt.show(vedo_sphere, vedo_fs, camera=cam)

elevation_amplitude = np.abs(diffraction_elevation_at_faces).max()
pressure_amplitude = np.abs(pressure_on_sphere).max()


def update(t):
    # Update the free surface elevation
    wave = np.real(diffraction_elevation_at_faces * np.exp(-1j * omega * t)).ravel()
    new_pts = vedo_fs.vertices
    new_pts[:, 2] = wave
    vedo_fs.vertices = new_pts
    vedo_fs.cmap("Blues_r", wave, vmin=-elevation_amplitude, vmax=elevation_amplitude)
    # Blues_r is used to color the crests in white and the troughs in blue

    # Update the pressure on the sphere
    pressure = np.real(pressure_on_sphere * np.exp(-1j * omega * t)).ravel()
    vedo_sphere.cmap(
        "RdBu_r", pressure, vmin=-pressure_amplitude, vmax=pressure_amplitude
    )
    vedo_sphere.add_scalarbar(title="dynamic pressure")


video_duration = 3 * T
video_fps = 24
video = Video("sphere_diffraction.mp4", duration=video_duration, fps=video_fps)
t_range = np.linspace(0.0, video_duration, int(video_duration * video_fps))
for t in t_range:
    update(t)
    plt.render()
    # video.add_frame()

video.close()
plt.close()
