from time import sleep
import numpy as np
import capytaine as cpt
from vedo import Plotter, Text2D, Mesh, Video
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation, airy_waves_pressure

omega = 3.0
T = 2*np.pi/omega
mesh = cpt.mesh_sphere(radius=3, center=(0, 0, 0), resolution=(20, 20))
sphere = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=["Heave"])).immersed_part()
solver = cpt.BEMSolver()
diffraction_problem = cpt.DiffractionProblem(body=sphere, wave_direction=0.0, omega=omega)
diffraction_result = solver.solve(diffraction_problem)

free_surface_mesh = cpt.mesh_rectangle(size=(50.0, 50.0), center=(0.0, 0.0, 0.0), resolution=(100, 100))

amplitude = 0.2
diffraction_elevation_at_faces = amplitude*(
        solver.compute_free_surface_elevation(free_surface_mesh.vertices, diffraction_result) +
        airy_waves_free_surface_elevation(free_surface_mesh.vertices, diffraction_problem)
        )
pressure_on_sphere = amplitude*(
        solver.compute_pressure(sphere.mesh.vertices, diffraction_result) +
        airy_waves_pressure(sphere.mesh.vertices, diffraction_problem)
        )


duration = 3*T
fps = 24
video = Video("vedo_video.mp4", duration=duration, fps=fps)

txt = Text2D(font='Brachium', pos='bottom-left', bg='yellow5')
vedo_sphere = Mesh([sphere.mesh.vertices, sphere.mesh.faces])
vedo_fs = Mesh([free_surface_mesh.vertices, free_surface_mesh.faces])

cam = dict(
    pos=(-30.0, -30.0, -30.0),
    viewup=(0, 0, 1),
)
plt = Plotter(axes=1, size=(1000,700), interactive=False)
plt.show(vedo_sphere, vedo_fs, txt, "diffracted wave", camera=cam)

t_range = np.linspace(0.0, duration, int(duration*fps))
for t in t_range:
    wave = np.real(diffraction_elevation_at_faces * np.exp(-1j*omega * t)).ravel()
    pressure = np.real(pressure_on_sphere * np.exp(-1j*omega * t)).ravel()
    vedo_fs.cmap("Blues", wave, vmin=-np.abs(diffraction_elevation_at_faces).max(), vmax=np.abs(diffraction_elevation_at_faces).max())
    vedo_sphere.cmap("RdBu_r", pressure, vmin=-np.abs(pressure_on_sphere).max(), vmax=np.abs(pressure_on_sphere).max())
    txt.text(f"time: {t/T:.2f} period")
    # newpts = grid.points
    # newpts[:, 2] = wave
    # grid.points = newpts
    new_pts = vedo_fs.vertices
    new_pts[:, 2] = wave
    vedo_fs.vertices = new_pts
    plt.render()
    video.add_frame()

video.close()
plt.close()
