#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

import logging

import numpy as np
from numpy.linalg import norm

from capytaine.reference_bodies import FreeSurface
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

        if isinstance(S, BlockToeplitzMatrix):
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

    # def get_potential_at_point(self, problem, point):
    #     if not problem.sources:
    #         raise Exception("Need to solve the problem")

    #     s = np.zeros((len(problem.body.faces_centers)), dtype=np.complex64)

    #     for i in range(problem.body.nb_faces):
    #         s[i] += -_Green.green_1.compute_s0(point,
    #                                            problem.body.vertices[problem.body.faces[i, :], :],
    #                                            problem.body.faces_centers[i, :],
    #                                            problem.body.faces_normals[i, :],
    #                                            problem.body.faces_areas[i],
    #                                            problem.body.faces_radiuses[i],
    #                                            )[0]/(4*np.pi)

    #         if problem.depth == np.infty:
    #             p_mirror = point.copy()
    #             p_mirror[2] = -p_mirror[2]
    #             s[i] += -_Green.green_1.compute_s0(p_mirror,
    #                                                problem.body.vertices[problem.body.faces[i, :], :],
    #                                                problem.body.faces_centers[i, :],
    #                                                problem.body.faces_normals[i, :],
    #                                                problem.body.faces_areas[i],
    #                                                problem.body.faces_radiuses[i],
    #                                                )[0]/(4*np.pi)

    #             s[i] += _Green.green_2.vnsinfd(problem.wavenumber, problem.body.faces_centers[i, :], point)[0]

    #         else:
    #             p_mirror = point.copy()
    #             p_mirror[2] = -2*problem.depth - p_mirror[2]
    #             s[i] += _Green.green_1.compute_s0(p_mirror,
    #                                                problem.body.vertices[problem.body.faces[i, :], :],
    #                                                problem.body.faces_centers[i, :],
    #                                                problem.body.faces_normals[i, :],
    #                                                problem.body.faces_areas[i],
    #                                                problem.body.faces_radiuses[i],
    #                                                )[0]/(4*np.pi)

    #             s[i] += _Green.green_2.vnsfd(problem.wavenumber, problem.body.faces_centers[i, :], point, problem.depth)[0]

    #     return s @ problem.sources["Heave"]

    def get_free_surface(self, problem):
        S, _ = FreeSurface().build_matrices(
            problem.body,
            free_surface=problem.free_surface,
            sea_bottom=problem.sea_bottom,
            wavenumber=problem.wavenumber
        )
        phi = S @ problem.sources["Heave"]
        return 1j*problem.omega/problem.g * phi.reshape(15, 15)

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)
