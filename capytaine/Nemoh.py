#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
# from itertools import chain, accumulate

import numpy as np

from capytaine.Toeplitz_matrices import (identity_matrix_of_same_shape_as, solve,
                                         BlockToeplitzMatrix, BlockCirculantMatrix)
from capytaine.symmetries import ReflectionSymmetry, TranslationalSymmetry, AxialSymmetry
from capytaine.tools.max_length_dict import MaxLengthDict
from capytaine.tools.exponential_decomposition import exponential_decomposition, error_exponential_decomposition
import capytaine._Green as _Green


LOG = logging.getLogger(__name__)


class Nemoh:
    """Solver for the BEM problem based on Nemoh's Green function."""

    def __init__(self, npinte=251, max_stored_exponential_decompositions=50):
        self.XR, self.XZ, self.APD = _Green.initialize_green_2.initialize_green(328, 46, npinte)
        LOG.info("Initialize Nemoh's Green function.")

        self.exponential_decompositions = MaxLengthDict(max_length=max_stored_exponential_decompositions)

        self.__cache__ = {'Green0': {}, 'Green1': {}, 'Green2': {}}

    def solve(self, problem, keep_details=False):
        """Solve the BEM problem using Nemoh."""

        LOG.info("Solve %s.", problem)

        if problem.depth < np.infty:
            self.compute_exponential_decomposition(problem)

        if problem.wavelength < 8*problem.body.mesh.faces_radiuses.max():
            LOG.warning(f"Resolution of the mesh (max_radius={problem.body.mesh.faces_radiuses.max():.2e}) "
                        f"might be insufficient for this wavelength (wavelength/8={problem.wavelength/8:.2e})!")

        S, V = self.build_matrices(
            problem.body, problem.body,
            free_surface=problem.free_surface, sea_bottom=problem.sea_bottom, wavenumber=problem.wavenumber
        )

        identity = identity_matrix_of_same_shape_as(V)
        sources = solve(V + identity/2, problem.boundary_condition)
        potential = S @ sources

        result = problem.make_results_container()
        if keep_details:
            result.sources = sources
            result.potential = potential

        for influenced_dof_name, influenced_dof in problem.body.dofs.items():
            integrated_potential = - problem.rho * potential @ (influenced_dof * problem.body.mesh.faces_areas)
            result.store_force(influenced_dof_name, integrated_potential)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and damping radiation coefficients.

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        pool = Pool(processes=processes)
        return pool.map(self.solve, problems)

    ####################
    #  Initialization  #
    ####################

    def compute_exponential_decomposition(self, pb):
        """Return the decomposition a part of the finite depth Green function as a sum of
        exponential functions.
        """

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
            self.exponential_decompositions.move_to_end(
                key=(pb.dimensionless_omega, pb.dimensionless_wavenumber), last=True)

    #######################
    #  Building matrices  #
    #######################

    def build_matrices(self, body1, body2,
                       free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0,
                       force_full_computation=False, _rec_depth=1):
        """Assemble the influence matrices.
        The method is basically an ugly multiple dispatch on the kind of bodies.
        For symmetric structures, the method is called recursively on the sub-bodies.

        Parameters
        ----------
        body1: FloatingBody
            receiving body (where the potential is measured)
        body2: FloatingBody
            source body (over which the source distribution is integrated)
        free_surface: float
            position of the free surface (default: z = 0)
        sea_bottom: float
            position of the sea bottom (default: z = -∞)
        wavenumber: float
            wavenumber (default: 1)
        force_full_computation: bool
            if True, the symmetries are NOT used to speed up the computation (default: False)
        _rec_depth: int
            internal parameter: recursion depth for pretty log printing

        Returns
        -------
        S: array of shape (body1.nb_faces, body2.nb_faces)
            influence matrix (integral of the Green function)
        V: array of shape (body1.nb_faces, body2.nb_faces)
            influence matrix (integral of the derivative of the Green function)
        """

        if (isinstance(body1, ReflectionSymmetry)
                and isinstance(body2, ReflectionSymmetry)
                and body1.plane == body2.plane
                and not force_full_computation):

            LOG.debug("\t"*_rec_depth +
                      f"Evaluating matrix of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"using mirror symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_a, V_a = self.build_matrices(body1.subbodies[0], body2.subbodies[0],
                                           free_surface, sea_bottom, wavenumber, force_full_computation, _rec_depth + 1)
            S_b, V_b = self.build_matrices(body1.subbodies[0], body2.subbodies[1],
                                           free_surface, sea_bottom, wavenumber, force_full_computation, _rec_depth + 1)

            return BlockToeplitzMatrix([S_a, S_b]), BlockToeplitzMatrix([V_a, V_b])

        elif (isinstance(body1, TranslationalSymmetry)
              and isinstance(body2, TranslationalSymmetry)
              and np.allclose(body1.translation, body2.translation)
              and body1.nb_subbodies == body2.nb_subbodies
              and not force_full_computation):

            LOG.debug("\t"*_rec_depth +
                      f"Evaluating matrix of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"using translational symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_list, V_list = [], []
            for subbody in body2.subbodies:
                S, V = self.build_matrices(body1.subbodies[0], subbody,
                                           free_surface, sea_bottom, wavenumber, force_full_computation, _rec_depth + 1)
                S_list.append(S)
                V_list.append(V)
            return BlockToeplitzMatrix(S_list), BlockToeplitzMatrix(V_list)

        elif (isinstance(body1, AxialSymmetry)
              and body1 is body2  # TODO: Generalize: if body1.axis == body2.axis
              and not force_full_computation):

            LOG.debug("\t" * _rec_depth +
                      f"Evaluating matrix of {body1.name} on itself "
                      f"using rotation symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_list, V_list = [], []
            for subbody in body2.subbodies[:body2.nb_subbodies // 2 + 1]:
                S, V = self.build_matrices(body1.subbodies[0], subbody,
                                           free_surface, sea_bottom, wavenumber, force_full_computation, _rec_depth + 1)
                S_list.append(S)
                V_list.append(V)

            if body1.nb_subbodies % 2 == 0:
                return BlockCirculantMatrix(S_list, size=body1.nb_subbodies), BlockCirculantMatrix(V_list, size=body1.nb_subbodies)
            else:
                return BlockCirculantMatrix(S_list, size=body1.nb_subbodies), BlockCirculantMatrix(V_list, size=body1.nb_subbodies)

        #   elif (isinstance(body1, CollectionOfFloatingBodies)):
        #     S = np.empty((body1.nb_faces, body2.nb_faces), dtype=np.complex64)
        #     V = np.empty((body1.nb_faces, body2.nb_faces), dtype=np.complex64)
        #
        #     nb_faces = list(accumulate(chain([0], (body.nb_faces for body in body1.subbodies))))
        #     for (i, j), body in zip(zip(nb_faces, nb_faces[1:]), body1.subbodies):
        #         matrix_slice = (slice(i, j), slice(None, None))
        #         S[matrix_slice], V[matrix_slice] = self.build_matrices(body1, body2, **kwargs)
        #
        #     return S, V

        else:
            LOG.debug("\t" * _rec_depth +
                      f"Evaluating matrix of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S = np.zeros((body1.mesh.nb_faces, body2.mesh.nb_faces), dtype=np.complex64)
            V = np.zeros((body1.mesh.nb_faces, body2.mesh.nb_faces), dtype=np.complex64)

            S0, V0 = self._build_matrices_0(body1, body2, _rec_depth)
            S += S0
            V += V0

            if free_surface < np.infty:

                S1, V1 = self._build_matrices_1(body1, body2, free_surface, sea_bottom, _rec_depth)
                S += S1
                V += V1

                S2, V2 = self._build_matrices_2(body1, body2, free_surface, sea_bottom, wavenumber, _rec_depth)
                S += S2
                V += V2

            return S, V

    def _build_matrices_0(self, body1, body2, _rec_depth=1):
        """Compute the first part of the influence matrices of self on body."""
        if body1 not in self.__cache__['Green0']:
            self.__cache__['Green0'][body1] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug("\t"*_rec_depth +
                      f"\tCreate Green0 cache (max_length={body1.nb_matrices_to_keep}) for {body1.name}")

        if body2 not in self.__cache__['Green0'][body1]:
            LOG.debug("\t"*_rec_depth +
                      f"\tComputing matrix 0 of {body1.name} on {'itself' if body2 is body1 else body2.name}")
            S0, V0 = _Green.green_1.build_matrix_0(
                body1.mesh.faces_centers, body1.mesh.faces_normals,
                body2.mesh.vertices,      body2.mesh.faces + 1,
                body2.mesh.faces_centers, body2.mesh.faces_normals,
                body2.mesh.faces_areas,   body2.mesh.faces_radiuses,
                )

            self.__cache__['Green0'][body1][body2] = (S0, V0)
        else:
            LOG.debug("\t"*_rec_depth +
                      f"\tRetrieving stored matrix 0 of {body1.name} on {'itself' if body2 is body1 else body2.name}")
            S0, V0 = self.__cache__['Green0'][body1][body2]

        return S0, V0

    def _build_matrices_1(self, body1, body2, free_surface, sea_bottom, _rec_depth=1):
        """Compute the second part of the influence matrices of body1 on body2."""
        if body1 not in self.__cache__['Green1']:
            self.__cache__['Green1'][body1] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug("\t"*_rec_depth +
                      f"\tCreate Green1 cache (max_length={body1.nb_matrices_to_keep}) for {body1.name}")

        depth = free_surface - sea_bottom
        if (body2, depth) not in self.__cache__['Green1'][body1]:
            LOG.debug("\t"*_rec_depth +
                      f"\tComputing matrix 1 of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"for depth={depth:.2e}")

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
                reflect_point(body1.mesh.faces_centers), reflect_vector(body1.mesh.faces_normals),
                body2.mesh.vertices,      body2.mesh.faces + 1,
                body2.mesh.faces_centers, body2.mesh.faces_normals,
                body2.mesh.faces_areas,   body2.mesh.faces_radiuses,
                )

            if depth == np.infty:
                self.__cache__['Green1'][body1][(body2, np.infty)] = (-S1, -V1)
                return -S1, -V1
            else:
                self.__cache__['Green1'][body1][(body2, depth)] = (S1, V1)
                return S1, V1
        else:
            S1, V1 = self.__cache__['Green1'][body1][(body2, depth)]
            LOG.debug("\t"*_rec_depth +
                      f"\tRetrieving stored matrix 1 of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"for depth={depth:.2e}")
            return S1, V1

    def _build_matrices_2(self, body1, body2, free_surface, sea_bottom, wavenumber, _rec_depth=1):
        """Compute the third part of the influence matrices of body1 on body2."""
        if body1 not in self.__cache__['Green2']:
            self.__cache__['Green2'][body1] = MaxLengthDict({}, max_length=body1.nb_matrices_to_keep)
            LOG.debug("\t"*_rec_depth +
                      f"\tCreate Green2 cache (max_length={body1.nb_matrices_to_keep}) for {body1.name}")

        depth = free_surface - sea_bottom
        if (body2, depth, wavenumber) not in self.__cache__['Green2'][body1]:
            LOG.debug("\t"*_rec_depth +
                      f"\tComputing matrix 2 of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"for depth={depth:.2e} and k={wavenumber:.2e}")
            if depth == np.infty:
                lamda_exp = np.empty(31, dtype=np.float32)
                a_exp = np.empty(31, dtype=np.float32)
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.mesh.faces_centers, body1.mesh.faces_normals,
                    body2.mesh.faces_centers, body2.mesh.faces_areas,
                    wavenumber,         0.0,
                    self.XR, self.XZ, self.APD,
                    lamda_exp, a_exp, n_exp,
                    body1 is body2
                    )
            else:
                # Get the last computed exponential decomposition.
                a_exp, lamda_exp = next(reversed(self.exponential_decompositions.values()))
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    body1.mesh.faces_centers, body1.mesh.faces_normals,
                    body2.mesh.faces_centers, body2.mesh.faces_areas,
                    wavenumber, depth,
                    self.XR, self.XZ, self.APD,
                    lamda_exp, a_exp, n_exp,
                    body1 is body2
                    )

            self.__cache__['Green2'][body1][(body2, depth, wavenumber)] = (S2, V2)
        else:
            S2, V2 = self.__cache__['Green2'][body1][(body2, depth, wavenumber)]
            LOG.debug("\t"*_rec_depth +
                      f"\tRetrieving stored matrix 2 of {body1.name} on {'itself' if body2 is body1 else body2.name} "
                      f"for depth={depth:.2e} and k={wavenumber:.2e}")

        return S2, V2

    #######################
    #  Compute potential  #
    #######################

    def get_potential_on_mesh(self, result, mesh):
        """Compute the potential on a mesh for the potential field of a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        mesh : FloatingBody
            a meshed floating body

        Returns
        -------
        array
            potential on the faces of the mesh
        """
        LOG.info(f"Compute potential on {mesh.name} for {result}.")

        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        S, _ = self.build_matrices(
            mesh,
            result.body,
            free_surface=result.free_surface,
            sea_bottom=result.sea_bottom,
            wavenumber=result.wavenumber
        )

        phi = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        free_surface : FloatingBody
            a meshed free surface

        Returns
        -------
        array
            the free surface elevation on each faces of the meshed free surface
        """
        return 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface)

