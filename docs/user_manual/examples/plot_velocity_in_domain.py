import numpy as np
import capytaine as cpt

# mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
mesh = cpt.mesh_horizontal_cylinder(radius=1.0, length=1.0, center=(0, 0, 0), resolution=(5, 20, 5)).immersed_part()
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
velocities = {}
solver = cpt.BEMSolver()
for dof in body.dofs:
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof=dof)
    res = solver.solve(pb)
    velocities[dof] = solver.get_velocity_on_mesh(res, mesh)

Y, Z = np.meshgrid(np.linspace(-3.0, 3.0, 40), np.linspace(-3.0, -0.05, 20))
points = np.vstack([np.zeros(Y.size), Y.ravel(), Z.ravel()]).T
points = np.array([[x, y, z] for (x, y, z) in points if y**2 + z**2 > 1])
velocities = {}
for dof in body.dofs:
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof=dof)
    res = solver.solve(pb)
    velocities[dof] = solver.get_velocity_at_points(res, points)

import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(12, 4))
for ax, dof in zip(axs.ravel(), ["Sway", "Heave"]):
    y_velocities = -velocities[dof][:, 1].real
    z_velocities = -velocities[dof][:, 2].real
    ax.quiver(points[:, 1], points[:, 2], y_velocities, z_velocities, angles="xy")
    ax.axis('equal')
    theta = np.linspace(-np.pi, 0.0, 100)
    ax.plot(np.cos(theta), np.sin(theta), color="red")
    ax.set_title(dof)
    ax.set_xlabel("y")
    ax.set_ylabel("z")
plt.show()


# import matplotlib.pyplot as plt
# fig, axs = plt.subplots(2, 3, subplot_kw=dict(projection="3d"))
# for ax, dof in zip(axs.ravel(), velocities):
#     ax.quiver(*zip(*mesh.faces_centers), *zip(*velocities[dof]), length=0.5)
#     ax.set_title(dof)
# plt.show()

