#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

import numpy as np
from numpy.linalg import solve

import capytaine._Green as _Green


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """
    def __init__(self):
        _Green.initialize_green_2.initialize_green()

    def build_matrices(self, pb):

        # for body1, body2 in product(bodies, repeat=2):
        body1 = pb.bodies[0]
        body2 = pb.bodies[0]

        S = np.zeros((body1.nb_faces, body2.nb_faces), dtype=np.complex64)
        V = np.zeros((body1.nb_faces, body2.nb_faces), dtype=np.complex64)

        S0, V0 = _Green.green_1.build_matrix_0(
            body1.faces_centers, body1.faces_normals,
            body2.vertices,      body2.faces + 1,
            body2.faces_centers, body2.faces_normals,
            body2.faces_areas,   body2.faces_radiuses,
            )
        S += S0
        V += V0

        if pb.free_surface < np.infty:

            if pb.sea_bottom == -np.infty:

                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*pb.free_surface - x[:, 2] 
                    return y
                def reflect_vector(x):
                    y = x.copy()
                    y[:, 2] = -x[:, 2]
                    return y

                S1, V1 = _Green.green_1.build_matrix_0(
                    reflect_point(body1.faces_centers), reflect_vector(body1.faces_normals),
                    body2.vertices,      body2.faces + 1,
                    body2.faces_centers, body2.faces_normals,
                    body2.faces_areas,   body2.faces_radiuses,
                    )

                S += -S1
                V += -V1

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,  
                    pb.wavenumber,       0.0
                    )

                S += S2
                V += V2
                
            else:

                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*pb.sea_bottom - x[:, 2]
                    return y
                def reflect_vector(x):
                    y = x.copy()
                    y[:, 2] = -x[:, 2]
                    return y

                S1, V1 = _Green.green_1.build_matrix_0(
                    reflect_point(body1.faces_centers), reflect_vector(body1.faces_normals),
                    body2.vertices,      body2.faces + 1,
                    body2.faces_centers, body2.faces_normals,
                    body2.faces_areas,   body2.faces_radiuses,
                    )

                S += S1
                V += V1

                _Green.initialize_green_2.lisc(
                    pb.omega**2*pb.depth/pb.g,
                    pb.wavenumber*pb.depth
                )

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,  
                    pb.wavenumber,       pb.depth
                    )

                S += S2
                V += V2

        return S, V

    def solve(self, problem):
        """Solve the BEM problem using Nemoh.
        Return the added mass and added damping.
        """

        S, V = self.build_matrices(problem)

        identity = np.identity(V.shape[0], dtype=np.float32)

        for body in problem.bodies:
            for dof in body.dof:
                sources = solve(V + identity/2, body.dof[dof])
                potential = S @ sources

                complex_coef = - problem.rho * potential @ \
                    (body.dof[dof] * body.faces_areas)

                added_mass = complex_coef.real
                added_damping = problem.omega * complex_coef.imag

        return added_mass, added_damping

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)
