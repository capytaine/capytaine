#!/usr/bin/env python

import numpy as np
import capytaine as cpt

# Initialize floating body
sphere = cpt.Sphere(
    radius=1.0,          # Dimension
    center=(0, 0, -2),   # Position
    nphi=20, ntheta=20,  # Fineness of the mesh
)

# DEFINE THE DOFS

# Manually defined heave,
# that is a vertical unit vector for all faces.
sphere.dofs["Heave"] = np.array(
    [(0, 0, 1) for face in sphere.mesh.faces]
)

# Bulging of the sphere,
# that is a deformation vector normal to the body at all faces.
# We can simply point to the normal vectors of the mesh for that.
sphere.dofs["Bulge"] = sphere.mesh.faces_normals

# Shearing of the sphere in the x direction.
# The deformation vector on a face is computed from the position of the face.
sphere.dofs["x-shear"] = np.array(
    [(np.cos(np.pi*z/2), 0, 0) for x, y, z in sphere.mesh.faces_centers]
)

# SOLVE DIFFRACTION PROBLEMS
solver = cpt.BEMSolver()

# Solve the problem for β=0 (incoming wave in the x direction).
problem_1 = cpt.DiffractionProblem(body=sphere, wave_direction=0, omega=1.0)
result_1 = solver.solve(problem_1)

# Solve the problem for β=π/2 (incoming wave in the y direction).
problem_2 = cpt.DiffractionProblem(body=sphere, wave_direction=np.pi/2, omega=1.0)
result_2 = solver.solve(problem_2)

# Print the generalized diffraction forces
# for the three dofs we defined
# for both values of the wave_direction β.
for result in [result_1, result_2]:
    print(f"Angle: {result.wave_direction:.2f}")
    for dof in sphere.dofs:
        force = result.forces[dof]
        print(f"{dof}: {np.abs(force):.2f}·exp({np.angle(force):.2f}i) N")
    print()

