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
from capytaine.tools.max_length_dict import MaxLengthDict
from capytaine.tools.exponential_decomposition import exponential_decomposition, error_exponential_decomposition


LOG = logging.getLogger(__name__)


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """

    def __init__(self):
        self.XR = _Green.initialize_green_2.initialize_green()
        LOG.info("Initialize Nemoh's Green function.")

        self.exponential_decompositions = MaxLengthDict(max_length=10)

    def solve(self, problem, keep_details=False):
        """Solve the BEM problem using Nemoh.
        """

        LOG.info("Solve %s.", problem)

        if problem.depth < np.infty:
            self.compute_exponential_decomposition(problem)

        S, V = problem.body.build_matrices(
            self,
            problem.body,
            free_surface=problem.free_surface,
            sea_bottom=problem.sea_bottom,
            wavenumber=problem.wavenumber
        )

        if isinstance(S, BlockCirculantMatrix):
            identity = block_circulant_identity(V.nb_blocks, V.block_size, dtype=np.float32)
        elif isinstance(S, BlockToeplitzMatrix):
            identity = block_Toeplitz_identity(V.nb_blocks, V.block_size, dtype=np.float32)
        else:
            identity = np.identity(V.shape[0], dtype=np.float32)

        sources = solve(V + identity/2, problem.boundary_condition)
        potential = S @ sources

        results = problem.make_results_container()

        if keep_details:
            results.keep_details = True
            results.S = S
            results.V = V
            results.sources = sources
            results.potential = potential

        for influenced_dof_name, influenced_dof in problem.body.dofs.items():
            force = - problem.rho * potential @ (influenced_dof * problem.body.faces_areas)
            results.store_force(influenced_dof_name, force)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and damping radiation coefficients.

        return results

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)

    ####################
    #  Initialization  #
    ####################
    def compute_exponential_decomposition(self, pb):
        """Return the decomposition a part of the finite depth Green function as a sum of
        exponentials."""

        LOG.debug(f"Initialize Nemoh's finite depth Green function for omega=%.2e and depth=%.2e", pb.omega, pb.depth)
        if (pb.dimensionless_omega, pb.dimensionless_wavenumber) not in self.exponential_decompositions:

            # The function that will be approximated.
            @np.vectorize
            def f(x):
                return _Green.initialize_green_2.ff(x, pb.dimensionless_omega,
                                                    pb.dimensionless_wavenumber)

            # Try different increasing number of exponentials
            for n_exp in range(4, 31, 2):

                # The coefficients are computed on a resolution of 4*n_exp+1 ...
                X = np.linspace(-0.1, 20.0, 4*n_exp+1)
                a, lamda = exponential_decomposition(X, f(X), n_exp)

                # ... and they are evaluated on a finer discretization.
                X = np.linspace(-0.1, 20.0, 8*n_exp+1)
                if error_exponential_decomposition(X, f(X), a, lamda) < 1e-4:
                    break

            else:
                LOG.warning(f"No suitable exponential decomposition has been found for {pb}.")

            # Convert to precision wanted by Fortran code.
            a = a.astype(np.float32)
            lamda = lamda.astype(np.float32)

            # Temporary trick: expand arrays to fix size hard-coded in Fortran module.
            a = np.r_[a, np.zeros(31-len(a), dtype=np.float32)]
            lamda = np.r_[lamda, np.zeros(31-len(lamda), dtype=np.float32)]

            self.exponential_decompositions[(pb.dimensionless_omega, pb.dimensionless_wavenumber)] = (a, lamda)

        else:
            self.exponential_decompositions.move_to_end(key=(pb.dimensionless_omega, pb.dimensionless_wavenumber), last=True)

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
                lamda_exp = np.empty(31, dtype=np.float32)
                a_exp = np.empty(31, dtype=np.float32)
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,
                    wavenumber,         0.0,
                    self.XR, lamda_exp, a_exp, n_exp,
                    body1 is body2
                    )
            else:
                # Get the last computed exponential decomposition.
                a_exp, lamda_exp= next(reversed(self.exponential_decompositions.values()))
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.faces_centers, body1.faces_normals,
                    body2.faces_centers, body2.faces_areas,
                    wavenumber,         depth,
                    self.XR, lamda_exp, a_exp, n_exp,
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

    def get_potential_on_mesh(self, solved_problem, mesh, dof=None):
        LOG.info(f"Compute potential on {mesh.name} for {solved_problem}.")

        if len(solved_problem.sources) == 0:
            if not solved_problem.keep_details:
                raise Exception(f"""The detail of the sources of {solved_problem} are not stored by the solver
                so they can't be used for later computation of the potential.
                Please run the solver with keep_details=True.""")
            else:
                raise Exception(f"{solved_problem} need to be solved with Nemoh.solve before computing potential anywhere.")

        S, _ = mesh.build_matrices(
            self,
            solved_problem.body,
            free_surface=solved_problem.free_surface,
            sea_bottom=solved_problem.sea_bottom,
            wavenumber=solved_problem.wavenumber
        )

        phi = S @ solved_problem.sources

        LOG.info(f"Done computing potential on {mesh.name} for {solved_problem}.")

        return phi

    def get_free_surface(self, solved_problem, free_surface, dof=None):
        return 1j*solved_problem.omega/solved_problem.g * self.get_potential_on_mesh(solved_problem, free_surface, dof=dof)

