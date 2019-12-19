#!/usr/bin/env python

import numpy as np
from numpy import dot, pi
from numpy.linalg import norm

import capytaine as cpt
from capytaine.green_functions.abstract_green_function import AbstractGreenFunction

class MyGreenFunction(AbstractGreenFunction):
    """An example of a custom routine to evaluate the Green function."""

    def evaluate(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
        """The main method that needs to be implemented in the class."""

        if free_surface == np.infty and sea_bottom == -np.infty:

            # Initialize the matrices
            S = np.zeros((mesh1.nb_faces, mesh2.nb_faces))
            K = np.zeros((mesh1.nb_faces, mesh2.nb_faces))

            for i in range(mesh1.nb_faces):
                for j in range(mesh2.nb_faces):
                    area = mesh1.faces_areas[i]
                    if i != j:
                        # The integral on the face is computed using the value at the center of the face.
                        r = mesh1.faces_centers[i, :] - mesh2.faces_centers[j, :]

                        S[j, i] = -area/(4*pi*norm(r))
                        K[j, i] = -area*dot(r, mesh1.faces_normals[j, :])/(4*pi*norm(r)**2)

                    else:
                        # The kernel is singular.
                        # When i == j, we can't use the approximation above for the integral.
                        # Below, we take a rough approximation.
                        S[j, i] = -area/(4*pi)*3
                        K[j, i] = -area/(4*pi)*3

            if mesh1 is mesh2:
                for i in range(mesh1.nb_faces):
                    K[i, i] += 1/2

            return S, K

        else:
            raise NotImplementedError

# Define a solver using the Green function above.
solver = cpt.BEMSolver(green_function=MyGreenFunction())

# Example of a problem
sphere = cpt.Sphere(
    radius=1.0,          # Dimension
    center=(0, 0, -2),   # Position
    nphi=4, ntheta=4,  # Fineness of the mesh
)
sphere.add_translation_dof(name="Heave")
problem = cpt.RadiationProblem(body=sphere, free_surface=np.infty, sea_bottom=-np.infty, radiating_dof='Heave')

result = solver.solve(problem)
print(result.added_masses)
