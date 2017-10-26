#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

import logging

import numpy as np

from capytaine.problems import RadiationProblem, DiffractionProblem
from capytaine.Toeplitz_matrices import BlockToeplitzMatrix, block_Toeplitz_identity, solve
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

        LOG.info("Solve %s.", problem)

        if problem.depth < np.infty:
            _Green.initialize_green_2.lisc(
                problem.omega**2*problem.depth/problem.g,
                problem.wavenumber*problem.depth
            )
            LOG.debug(f"Initialize Nemoh's finite depth Green function for omega=%.2e and depth=%.2e", problem.omega, problem.depth)

        added_masses, added_dampings = [], []

        S, V = problem.body.build_matrices(
            problem.body,
            free_surface=problem.free_surface,
            sea_bottom=problem.sea_bottom,
            wavenumber=problem.wavenumber
        )

        if isinstance(S, BlockToeplitzMatrix):
            identity = block_Toeplitz_identity(V.nb_blocks, V.block_size, dtype=np.float32)
        else:
            identity = np.identity(V.shape[0], dtype=np.float32)

        if isinstance(problem, RadiationProblem):
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

        elif isinstance(problem, DiffractionProblem):
            normal_velocities = -(problem.Airy_wave(problem.body.faces_centers) *
                                  problem.body.faces_normals
                                  ).sum(axis=1)
            sources = solve(V + identity/2, normal_velocities)
            potential = S @ sources

            forces = []
            for _, influenced_dof in problem.body.dofs.items():
                force = - problem.rho * \
                    potential @ (influenced_dof * problem.body.faces_areas)
                forces.append(force)

            return np.array(forces)

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)
