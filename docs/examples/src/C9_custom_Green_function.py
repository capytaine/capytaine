import numpy as np
from numpy import dot, pi
from numpy.linalg import norm

import capytaine as cpt
from capytaine.green_functions.abstract_green_function import AbstractGreenFunction


class MyGreenFunction(AbstractGreenFunction):
    """An example of a custom routine to evaluate the Green function."""

    exportable_settings = {"green_function": "MyGreenFunction"}
    # `exportable_settings` is meant to contain metadata to be stored in output
    # dataset for reproductibility.

    floating_point_precision = "float64"  # Optional, could allow for RAM usage optimisations

    def evaluate(
        self,
        mesh1,
        mesh2,
        free_surface,
        water_depth,
        wavenumber,
        adjoint_double_layer=True,
        early_dot_product=True,
    ):
        """The main method that needs to be implemented in the class."""

        if not adjoint_double_layer:
            raise NotImplementedError(
                """MyGreenFunction only implements the adjoint double layer"""
            )

        if free_surface != np.inf and water_depth != np.inf:
            raise NotImplementedError(
                """MyGreenFunction only implements the case without a free surface"""
            )

        # Initialize the matrices
        S = np.zeros((mesh1.nb_faces, mesh2.nb_faces))
        if early_dot_product:
            gradient_shape = (mesh1.nb_faces, mesh2.nb_faces)
        else:
            gradient_shape = (mesh1.nb_faces, mesh2.nb_faces, 3)
        K = np.zeros(gradient_shape)

        for i in range(mesh1.nb_faces):
            for j in range(mesh2.nb_faces):
                area = mesh1.faces_areas[j]
                if i == j:
                    # Approximation of the integral of the 1/r Green
                    # function over the panel by taking the integral over a
                    # circle of same center and same area.
                    S[i, i] = -1 / (4 * pi) * 2 * np.sqrt(area * pi)
                    if mesh1 is mesh2:
                        if early_dot_product:
                            K[i, i] = 1 / 2
                        else:
                            K[i, i, :] = 1 / 2 * mesh1.faces_normals[i, :]
                    else:
                        if early_dot_product:
                            K[i, i] = 0.0
                        else:
                            K[i, i, :] = 0.0

                else:
                    # The integral on the face is computed using the value
                    # at the center of the face.
                    r = mesh1.faces_centers[i, :] - mesh2.faces_centers[j, :]

                    S[i, j] = -area / (4 * pi * norm(r))
                    gradient_rankine = area * r / (4 * pi * norm(r) ** 3)
                    if early_dot_product:
                        K[i, j] = dot(gradient_rankine, mesh1.faces_normals[i, :])
                    else:
                        K[i, j, :] = gradient_rankine

        return S, K


# Define a solver using the Green function above.
custom_solver = cpt.BEMSolver(green_function=MyGreenFunction())
default_solver = cpt.BEMSolver()  # for comparison

# Example of a problem
radius = 1.4
body = cpt.FloatingBody(
    mesh=cpt.mesh_sphere(radius=radius, resolution=(20, 20)),
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
)
problem = cpt.RadiationProblem(
    body=body, free_surface=np.inf, water_depth=np.inf, radiating_dof="Surge"
)

custom_result = custom_solver.solve(problem)
default_result = default_solver.solve(problem)

print("Added mass:")
print(f"    Analytical value: {2/3*pi*problem.rho*radius**3:.2e}")
print(f"    Custom solver:    {custom_result.added_masses['Surge']:.2e}")
print(f"    Default solver:   {default_result.added_masses['Surge']:.2e}")
# The analytical solution is known for this particular test case.
# The custom solver might be more accurate. That is probably a coincidence.
# Both should converge to the same solution when the mesh is refined.

print()
print("Timer:")
print(f"    Custom solver:  {custom_solver.timer['  Green function'].mean:.1e} s")
print(f"    Default solver: {default_solver.timer['  Green function'].mean:.1e} s")
