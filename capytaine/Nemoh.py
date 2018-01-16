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
from capytaine.tools import MaxLengthDict


LOG = logging.getLogger(__name__)


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """

    def __init__(self):
        self.XR = _Green.initialize_green_2.initialize_green()
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
            self,
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


    #######################
    #  Building matrices  #
    #######################

    def _build_matrices_0(self, body1, body2):
        """Compute the first part of the influence matrices of self on body."""
        if 'Green0' not in body1.__internals__:
            body1.__internals__['Green0'] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug(f"\t\tCreate Green0 dict (max_length={body1.nb_matrices_to_keep}) in {body1.name}")

        if body2 not in body1.__internals__['Green0']:
            LOG.debug(f"\t\tComputing matrix 0 of {body1.name} on {body2.name}")
            S0, V0 = _Green.green_1.build_matrix_0(
                body1.faces_centers, body1.faces_normals,
                body2.vertices,      body2.faces + 1,
                body2.faces_centers, body2.faces_normals,
                body2.faces_areas,   body2.faces_radiuses,
                )

            body1.__internals__['Green0'][body2] = (S0, V0)
        else:
            LOG.debug(f"\t\tRetrieving stored matrix 0 of {body1.name} on {body2.name}")
            S0, V0 = body1.__internals__['Green0'][body2]

        return S0, V0

    def _build_matrices_1(self, body1, body2, free_surface, sea_bottom):
        """Compute the second part of the influence matrices of body1 on body2."""
        if 'Green1' not in body1.__internals__:
            body1.__internals__['Green1'] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug(f"\t\tCreate Green1 dict (max_length={body1.nb_matrices_to_keep}) in {body1.name}")

        depth = free_surface - sea_bottom
        if (body2, depth) not in body1.__internals__['Green1']:
            LOG.debug(f"\t\tComputing matrix 1 of {body1.name} on {body2.name} for depth={depth:.2e}")
            def reflect_vector(x):
                y = x.copy()
                y[:, 2] = -x[:, 2]
                return y

            if depth == np.infty:
                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*free_surface - x[:, 2]
                    return y
            else:
                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*sea_bottom - x[:, 2]
                    return y

            S1, V1 = _Green.green_1.build_matrix_0(
                reflect_point(body1.faces_centers), reflect_vector(body1.faces_normals),
                body2.vertices,      body2.faces + 1,
                body2.faces_centers, body2.faces_normals,
                body2.faces_areas,   body2.faces_radiuses,
                )

            if depth == np.infty:
                body1.__internals__['Green1'][(body2, np.infty)] = (-S1, -V1)
                return -S1, -V1
            else:
                body1.__internals__['Green1'][(body2, depth)] = (S1, V1)
                return S1, V1
        else:
            S1, V1 = body1.__internals__['Green1'][(body2, depth)]
            LOG.debug(f"\t\tRetrieving stored matrix 1 of {body1.name} on {body2.name} for depth={depth:.2e}")
            return S1, V1

    def _build_matrices_2(self, body1, body2, free_surface, sea_bottom, wavenumber):
        """Compute the third part of the influence matrices of body1 on body2."""
        if 'Green2' not in body1.__internals__:
            body1.__internals__['Green2'] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug(f"\t\tCreate Green2 dict (max_length={body1.nb_matrices_to_keep}) in {body1.name}")

        depth = free_surface - sea_bottom
        if (body2, depth, wavenumber) not in body1.__internals__['Green2']:
            LOG.debug(f"\t\tComputing matrix 2 of {body1.name} on {body2.name} for depth={depth:.2e} and k={wavenumber:.2e}")
            if depth == np.infty:
                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,
                    wavenumber,         0.0,
                    self.XR,
                    body1 is body2
                    )
            else:
                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,
                    wavenumber,         depth,
                    self.XR,
                    body1 is body2
                    )

            body1.__internals__['Green2'][(body2, depth, wavenumber)] = (S2, V2)
        else:
            S2, V2 = body1.__internals__['Green2'][(body2, depth, wavenumber)]
            LOG.debug(f"\t\tRetrieving stored matrix 2 of {body1.name} on {body2.name} for depth={depth:.2e} and k={wavenumber:.2e}")

        return S2, V2

    #######################
    #  Compute potential  #
    #######################

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
            self,
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

