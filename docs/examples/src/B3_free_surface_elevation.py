import numpy as np
import capytaine as cpt
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import matplotlib.pyplot as plt

mesh = cpt.mesh_sphere(radius=1, resolution=(10, 10)).immersed_part()
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())

pb = cpt.DiffractionProblem(body=body, wave_direction=np.pi/2, wavelength=3.0)
solver = cpt.BEMSolver()
res = solver.solve(pb)

grid = np.meshgrid(np.linspace(-10.0, 10.0, 100), np.linspace(-10.0, 10.0, 100))
fse = solver.compute_free_surface_elevation(grid, res)
incoming_fse = airy_waves_free_surface_elevation(grid, res)

fig, ax = plt.subplots(1, 1)
ax.pcolormesh(grid[0], grid[1], np.real(fse + incoming_fse))
ax.add_patch(plt.Circle((0, 0), radius=1, color="black"))
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
ax.set_title("Wave crests and troughs at t=0 seen from above")
plt.show()
