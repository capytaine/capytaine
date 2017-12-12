#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

import logging

import numpy as np

from capytaine.problems import RadiationProblem, DiffractionProblem
from capytaine.Toeplitz_matrices import (BlockCirculantMatrix, block_circulant_identity,
                                         BlockToeplitzMatrix, block_Toeplitz_identity,
                                         solve)
import capytaine._Green as _Green


LOG = logging.getLogger(__name__)


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """
    def __init__(self):
        _Green.initialize_green_2.initialize_green()
        LOG.info("Initialize Nemoh's Green function.")

    def solve(self, problem, keep_details=False):
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

        if keep_details:
            problem.S = S
            problem.V = V

        if isinstance(S, BlockCirculantMatrix):
            identity = block_circulant_identity(V.nb_blocks, V.block_size, dtype=np.float32)
        elif isinstance(S, BlockToeplitzMatrix):
            identity = block_Toeplitz_identity(V.nb_blocks, V.block_size, dtype=np.float32)
        else:
            identity = np.identity(V.shape[0], dtype=np.float32)

        if isinstance(problem, RadiationProblem):
            for dof_name, radiating_dof in problem.body.dofs.items():
                sources = solve(V + identity/2, radiating_dof)
                potential = S @ sources

                if keep_details:
                    problem.sources[dof_name] = sources
                    problem.potential[dof_name] = potential

                for _, influenced_dof in problem.body.dofs.items():
                    complex_coef = - problem.rho * \
                        potential @ (influenced_dof * problem.body.faces_areas)

                    added_masses.append(complex_coef.real)
                    added_dampings.append(problem.omega * complex_coef.imag)

            LOG.info("Problem solved!")

            return np.array(added_masses).reshape((problem.body.nb_dofs, problem.body.nb_dofs)), \
                   np.array(added_dampings).reshape((problem.body.nb_dofs, problem.body.nb_dofs))

        elif isinstance(problem, DiffractionProblem):
            normal_velocities = -(problem.Airy_wave_velocity(problem.body.faces_centers) *
                                  problem.body.faces_normals
                                  ).sum(axis=1)
            sources = solve(V + identity/2, normal_velocities)
            potential = S @ sources

            if keep_details:
                problem.sources = sources
                problem.potential = potential

            forces = []
            for _, influenced_dof in problem.body.dofs.items():
                force = - problem.rho * \
                    potential @ (influenced_dof * problem.body.faces_areas)
                forces.append(force)

            LOG.info("Problem solved!")

            return np.array(forces)

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)

    def get_potential_on_mesh(self, problem, mesh, dof=None):
        LOG.info(f"Compute potential on {mesh.name} for {problem}.")

        if len(problem.sources) == 0:
            if not problem.keep_details:
                raise Exception(f"""The detail of the sources of {problem} are not stored by the solver 
                so they can't be used for later computation of the potential. 
                Please run the solver with keep_details=True.""")
            else:
                raise Exception(f"{problem} need to be solved with Nemoh.solve before computing potential anywhere.")

        S, _ = mesh.build_matrices(
            problem.body,
            free_surface=problem.free_surface,
            sea_bottom=problem.sea_bottom,
            wavenumber=problem.wavenumber
        )

        if isinstance(problem, RadiationProblem):
            if dof is None:
                raise Exception("Please chose a degree of freedom.")
            else:
                phi = S @ problem.sources[dof]
        elif isinstance(problem, DiffractionProblem):
            phi = S @ problem.sources

        LOG.info(f"Done computing potential on {mesh.name} for {problem}.")

        return phi

    def get_free_surface(self, problem, free_surface, dof=None):
        return 1j*problem.omega/problem.g * self.get_potential_on_mesh(problem, free_surface, dof=dof)

