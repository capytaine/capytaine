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

    def build_matrices(self, bodies, wavenumber, omega, depth, g):
        if depth < np.infty:
            _Green.initialize_green_2.lisc(omega**2*depth/g, wavenumber*depth)
        else:
            depth = 0.0

        # for body1, body2 in product(bodies, repeat=2):
        body1 = bodies[0]
        body2 = bodies[0]

        return _Green.build_matrices(
            body1.faces_centers, body1.faces_normals,
            body2.vertices,      body2.faces + 1,
            body2.faces_centers, body2.faces_normals,
            body2.faces_areas,   body2.faces_radiuses,
            wavenumber,          depth
            )

    def solve(self, problem):
        """Solve the BEM problem using Nemoh.
        Return the added mass and added damping.
        """

        S, V = self.build_matrices(
            problem.bodies, problem.wavenumber,
            problem.omega,  problem.depth,
            problem.g
        )

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
