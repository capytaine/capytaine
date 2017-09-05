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

        added_masses, added_dampings = [], []

        for radiating_body in problem.bodies:
            S, V = radiating_body.build_matrices(
                radiating_body,
                problem.free_surface,
                problem.sea_bottom,
                problem.wavenumber
            )

            identity = np.identity(V.shape[0], dtype=np.float32)

            for dof_name, radiating_dof in radiating_body.dof.items():
                sources = solve(V + identity/2, radiating_dof)

                for influenced_body in problem.bodies:
                    S, V = radiating_body.build_matrices(
                        influenced_body,
                        problem.free_surface,
                        problem.sea_bottom,
                        problem.wavenumber
                    )
                    potential = S @ sources

                    for dof_name, influenced_dof in influenced_body.dof.items():

                        complex_coef = - problem.rho * \
                            potential @ (influenced_dof * influenced_body.faces_areas)

                        added_masses.append(complex_coef.real)
                        added_dampings.append(problem.omega * complex_coef.imag)

        return np.array(added_masses).reshape((problem.nb_total_dofs, problem.nb_total_dofs)), \
               np.array(added_dampings).reshape((problem.nb_total_dofs, problem.nb_total_dofs))

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)
