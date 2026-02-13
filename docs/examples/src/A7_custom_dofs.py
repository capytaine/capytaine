import numpy as np
import capytaine as cpt

mesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), faces_max_radius=0.2)

# Initialize floating body
body = cpt.FloatingBody(mesh, lid_mesh=mesh.generate_lid())

# DEFINE THE DOFS

# Manually defined heave,
# that is a vertical unit vector for all faces.
body.dofs["Heave"] = np.array([(0, 0, 1) for _ in body.mesh.faces_centers])

# Bulging of the body,
# that is a deformation vector normal to the body at all faces.
# We can simply point to the normal vectors of the mesh for that.
body.dofs["Bulge"] = body.mesh.faces_normals

# Shearing of the body in the x direction.
# The deformation vector on a face is computed from the position of the face.
body.dofs["x-shear"] = np.array(
    [(np.cos(np.pi * z / 2), 0, 0) for x, y, z in body.mesh.faces_centers]
)

# SOLVE DIFFRACTION PROBLEMS
solver = cpt.BEMSolver()

# Solve the problem for β=0 (incoming wave in the x direction).
problem_1 = cpt.DiffractionProblem(body=body, wave_direction=0, omega=1.0)
result_1 = solver.solve(problem_1)

# Solve the problem for β=π/2 (incoming wave in the y direction).
problem_2 = cpt.DiffractionProblem(body=body, wave_direction=np.pi / 2, omega=1.0)
result_2 = solver.solve(problem_2)

# Print the generalized diffraction forces
# for the three dofs we defined
# for both values of the wave_direction β.
for result in [result_1, result_2]:
    print(f"Angle: {result.wave_direction:.2f}")
    for dof in body.dofs:
        force = result.forces[dof]
        print(f"{dof}: {np.abs(force):.2f}·exp({np.angle(force):.2f}i) N")
    print()
