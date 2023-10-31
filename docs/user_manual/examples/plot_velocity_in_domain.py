import numpy as np
import capytaine as cpt

# Solve the potential flow problems
radius = 1.0
mesh = cpt.mesh_horizontal_cylinder(radius=radius, length=1.0, center=(0, 0, 0), resolution=(5, 20, 5)).immersed_part()
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
solver = cpt.BEMSolver()
results = solver.solve_all([cpt.RadiationProblem(body=body, omega=1.0, radiating_dof=dof) for dof in body.dofs])

# Define a set of points in the fluid on which the fluid velocities will be computed
Y, Z = np.meshgrid(np.linspace(-3.0, 3.0, 40), np.linspace(-3.0, -0.05, 20))
points = np.vstack([np.zeros(Y.size), Y.ravel(), Z.ravel()]).T
# Filter out points that are inside the floating body
points = np.array([[x, y, z] for (x, y, z) in points if y**2 + z**2 > radius**2])

# Compute the velocity field for each of the result object, that is each of the radiating dof
velocities = {}
for res in results:
    velocities[res.radiating_dof] = solver.compute_velocity( points,res)

import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(12, 4))
for ax, dof in zip(axs.ravel(), ["Sway", "Heave"]):
    # Display the arrows of the velocity
    y_velocities = -velocities[dof][:, 1].real
    z_velocities = -velocities[dof][:, 2].real
    ax.quiver(points[:, 1], points[:, 2], y_velocities, z_velocities, angles="xy")
    ax.axis('equal')

    # Display the boundary of the floating body
    theta = np.linspace(-np.pi, 0.0, 100)
    ax.plot(np.cos(theta), np.sin(theta), color="red")
    ax.set_title(dof)
    ax.set_xlabel("y")
    ax.set_ylabel("z")

plt.show()
