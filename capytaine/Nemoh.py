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

    def solve(self, problem):
        """Solve the BEM problem using Nemoh.
        Return the added mass and added damping.
        """

        if problem.depth < np.infty:
            _Green.initialize_green_2.lisc(
                problem.omega**2*problem.depth/problem.g,
                problem.wavenumber*problem.depth
            )

        body = problem.bodies[0]
        S, V = body.build_matrices(problem.free_surface, problem.sea_bottom, problem.wavenumber)

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
