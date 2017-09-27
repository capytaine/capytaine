#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

import logging

import numpy as np
from numpy.linalg import solve

import capytaine._Green as _Green


LOG = logging.getLogger(__name__)


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """
    def __init__(self):
        _Green.initialize_green_2.initialize_green()
        LOG.info("Initialize Nemoh's Green function.")

    def solve(self, problem):
        """Solve the BEM problem using Nemoh.
        Return the added mass and added damping.
        """

        LOG.info(f"Solve {problem}.")

        if problem.depth < np.infty:
            _Green.initialize_green_2.lisc(
                problem.omega**2*problem.depth/problem.g,
                problem.wavenumber*problem.depth
            )
            LOG.info("Initialize Nemoh's finite depth Green function for omega={omega:.2e} and depth={depth:.2e}.")

        added_masses, added_dampings = [], []

        S, V = problem.body.build_matrices(
            problem.body,
            free_surface=problem.free_surface,
            sea_bottom=problem.sea_bottom,
            wavenumber=problem.wavenumber
        )

        identity = np.identity(V.shape[0], dtype=np.float32)

        for _, radiating_dof in problem.body.dofs.items():
            sources = solve(V + identity/2, radiating_dof)
            potential = S @ sources

            for _, influenced_dof in problem.body.dofs.items():
                complex_coef = - problem.rho * \
                    potential @ (influenced_dof * problem.body.faces_areas)

                added_masses.append(complex_coef.real)
                added_dampings.append(problem.omega * complex_coef.imag)

        return np.array(added_masses).reshape((problem.body.nb_dofs, problem.body.nb_dofs)), \
               np.array(added_dampings).reshape((problem.body.nb_dofs, problem.body.nb_dofs))

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)
