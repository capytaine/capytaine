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

FLOAT_PRECISION = np.float64


class Nemoh:
    """Solver for the BEM problem based on Nemoh's Green function."""

    def __init__(self, keep_matrices=True, npinte=251, max_stored_exponential_decompositions=50):
        self.XR, self.XZ, self.APD = _Green.initialize_green_2.initialize_green(328, 46, npinte)
        LOG.info("Initialize Nemoh's Green function.")

        self.exponential_decompositions = MaxLengthDict(max_length=max_stored_exponential_decompositions)

        self.keep_matrices = keep_matrices
        if self.keep_matrices:
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
            problem.body.mesh, problem.body.mesh,
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
            influenced_dof = np.sum(influenced_dof * problem.body.mesh.faces_normals, axis=1)
            integrated_potential = - problem.rho * potential @ (influenced_dof * problem.body.mesh.faces_areas)
            result.store_force(influenced_dof_name, integrated_potential)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and damping radiation coefficients.

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, processes=1):
        from multiprocessing import Pool
        with Pool(processes=processes) as pool:
            results = pool.map(self.solve, problems)
        return results

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
            a = a.astype(FLOAT_PRECISION)
            lamda = lamda.astype(FLOAT_PRECISION)

            # Temporary trick: expand arrays to fixed size hard-coded in Fortran module.
            a = np.r_[a, np.zeros(31-len(a), dtype=FLOAT_PRECISION)]
            lamda = np.r_[lamda, np.zeros(31-len(lamda), dtype=FLOAT_PRECISION)]

            self.exponential_decompositions[(pb.dimensionless_omega, pb.dimensionless_wavenumber)] = (a, lamda)

        else:
            self.exponential_decompositions.move_to_end(
                key=(pb.dimensionless_omega, pb.dimensionless_wavenumber), last=True)

    #######################
    #  Building matrices  #
    #######################

    def build_matrices(self, mesh1, mesh2,
                       free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0,
                       force_full_computation=False, _rec_depth=(1,)):
        """Assemble the influence matrices.
        The method is basically an ugly multiple dispatch on the kind of bodies.
        For symmetric structures, the method is called recursively on the sub-bodies.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float
            position of the free surface (default: z = 0)
        sea_bottom: float
            position of the sea bottom (default: z = -∞)
        wavenumber: float
            wavenumber (default: 1)
        force_full_computation: bool
            if True, the symmetries are NOT used to speed up the computation (default: False)
        _rec_depth: tuple
            internal parameter: recursion accumulator for pretty log printing and cache sizing

        Returns
        -------
        S: array of shape (mesh1.nb_faces, body2.nb_faces)
            influence matrix (integral of the Green function)
        V: array of shape (mesh1.nb_faces, body2.nb_faces)
            influence matrix (integral of the derivative of the Green function)
        """

        if (isinstance(mesh1, ReflectionSymmetry)
                and isinstance(mesh2, ReflectionSymmetry)
                and mesh1.plane == mesh2.plane
                and not force_full_computation):

            LOG.debug("\t" * len(_rec_depth) +
                      f"Evaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"using mirror symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_a, V_a = self.build_matrices(mesh1.submeshes[0], mesh2.submeshes[0],
                                           free_surface, sea_bottom, wavenumber,
                                           force_full_computation, _rec_depth + (2,))
            S_b, V_b = self.build_matrices(mesh1.submeshes[0], mesh2.submeshes[1],
                                           free_surface, sea_bottom, wavenumber,
                                           force_full_computation, _rec_depth + (2,))

            return BlockToeplitzMatrix([S_a, S_b]), BlockToeplitzMatrix([V_a, V_b])

        elif (isinstance(mesh1, TranslationalSymmetry)
              and isinstance(mesh2, TranslationalSymmetry)
              and np.allclose(mesh1.translation, mesh2.translation)
              and mesh1.nb_submeshes == mesh2.nb_submeshes
              and not force_full_computation):

            LOG.debug("\t" * len(_rec_depth) +
                      f"Evaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"using translational symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_list, V_list = [], []
            for subbody in mesh2.submeshes:
                S, V = self.build_matrices(mesh1.submeshes[0], subbody,
                                           free_surface, sea_bottom, wavenumber,
                                           force_full_computation, _rec_depth + (mesh2.nb_submeshes,))
                S_list.append(S)
                V_list.append(V)
            return BlockToeplitzMatrix(S_list), BlockToeplitzMatrix(V_list)

        elif (isinstance(mesh1, AxialSymmetry)
              and mesh1 is mesh2  # TODO: Generalize: if mesh1.axis == mesh2.axis
              and not force_full_computation):

            LOG.debug("\t" * len(_rec_depth) +
                      f"Evaluating matrix of {mesh1.name} on itself "
                      f"using rotation symmetry "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S_list, V_list = [], []
            for subbody in mesh2.submeshes[:mesh2.nb_submeshes // 2 + 1]:
                S, V = self.build_matrices(mesh1.submeshes[0], subbody,
                                           free_surface, sea_bottom, wavenumber,
                                           force_full_computation, _rec_depth + (mesh2.nb_submeshes // 2 + 1,))
                S_list.append(S)
                V_list.append(V)

            if mesh1.nb_submeshes % 2 == 0:
                return BlockCirculantMatrix(S_list, size=mesh1.nb_submeshes), BlockCirculantMatrix(V_list, size=mesh1.nb_submeshes)
            else:
                return BlockCirculantMatrix(S_list, size=mesh1.nb_submeshes), BlockCirculantMatrix(V_list, size=mesh1.nb_submeshes)

        #   elif (isinstance(mesh1, CollectionOfMeshes)):
        #     S = np.empty((mesh1.nb_faces, mesh2.nb_faces), dtype=np.complex64)
        #     V = np.empty((mesh1.nb_faces, mesh2.nb_faces), dtype=np.complex64)
        #
        #     nb_faces = list(accumulate(chain([0], (body.nb_faces for body in mesh1.submeshes))))
        #     for (i, j), body in zip(zip(nb_faces, nb_faces[1:]), mesh1.submeshes):
        #         matrix_slice = (slice(i, j), slice(None, None))
        #         S[matrix_slice], V[matrix_slice] = self.build_matrices(mesh1, mesh2, **kwargs)
        #
        #     return S, V

        else:
            LOG.debug("\t" * len(_rec_depth) +
                      f"Evaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

            S = np.zeros((mesh1.nb_faces, mesh2.nb_faces), dtype=np.complex64)
            V = np.zeros((mesh1.nb_faces, mesh2.nb_faces), dtype=np.complex64)

            S0, V0 = self._build_matrices_0(mesh1, mesh2, _rec_depth)
            S += S0
            V += V0

            if free_surface < np.infty:

                S1, V1 = self._build_matrices_1(mesh1, mesh2, free_surface, sea_bottom, _rec_depth)
                S += S1
                V += V1

                S2, V2 = self._build_matrices_2(mesh1, mesh2, free_surface, sea_bottom, wavenumber, _rec_depth)
                S += S2
                V += V2

            return S, V

    def _build_matrices_0(self, mesh1, mesh2, _rec_depth=(1,)):
        """Compute the first part of the influence matrices of self on body."""
        if self.keep_matrices and mesh1 not in self.__cache__['Green0']:
            self.__cache__['Green0'][mesh1] = MaxLengthDict({}, max_length=int(np.product(_rec_depth)))
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tCreate Green0 cache (max_length={int(np.product(_rec_depth))}) for {mesh1.name}")

        if not self.keep_matrices or mesh2 not in self.__cache__['Green0'][mesh1]:
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tComputing matrix 0 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name}")
            S0, V0 = _Green.green_1.build_matrix_0(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.vertices,      mesh2.faces + 1,
                mesh2.faces_centers, mesh2.faces_normals,
                mesh2.faces_areas,   mesh2.faces_radiuses,
                )

            if self.keep_matrices:
                self.__cache__['Green0'][mesh1][mesh2] = (S0, V0)

        else:
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tRetrieving stored matrix 0 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name}")
            S0, V0 = self.__cache__['Green0'][mesh1][mesh2]

        return S0, V0

    def _build_matrices_1(self, mesh1, mesh2, free_surface, sea_bottom, _rec_depth=(1,)):
        """Compute the second part of the influence matrices of mesh1 on mesh2."""
        if self.keep_matrices and mesh1 not in self.__cache__['Green1']:
            self.__cache__['Green1'][mesh1] = MaxLengthDict({}, max_length=int(np.product(_rec_depth)))
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tCreate Green1 cache (max_length={int(np.product(_rec_depth))}) for {mesh1.name}")

        depth = free_surface - sea_bottom
        if not self.keep_matrices or (mesh2, depth) not in self.__cache__['Green1'][mesh1]:
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tComputing matrix 1 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
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
                reflect_point(mesh1.faces_centers), reflect_vector(mesh1.faces_normals),
                mesh2.vertices,      mesh2.faces + 1,
                mesh2.faces_centers, mesh2.faces_normals,
                mesh2.faces_areas,   mesh2.faces_radiuses,
                )

            if depth == np.infty:
                if self.keep_matrices:
                    self.__cache__['Green1'][mesh1][(mesh2, np.infty)] = (-S1, -V1)
                return -S1, -V1
            else:
                if self.keep_matrices:
                    self.__cache__['Green1'][mesh1][(mesh2, depth)] = (S1, V1)
                return S1, V1

        else:
            S1, V1 = self.__cache__['Green1'][mesh1][(mesh2, depth)]
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tRetrieving stored matrix 1 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"for depth={depth:.2e}")
            return S1, V1

    def _build_matrices_2(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber, _rec_depth=(1,)):
        """Compute the third part of the influence matrices of mesh1 on mesh2."""
        if self.keep_matrices and mesh1 not in self.__cache__['Green2']:
            self.__cache__['Green2'][mesh1] = MaxLengthDict({}, max_length=int(np.product(_rec_depth)))
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tCreate Green2 cache (max_length={int(np.product(_rec_depth))}) for {mesh1.name}")

        depth = free_surface - sea_bottom
        if not self.keep_matrices or (mesh2, depth, wavenumber) not in self.__cache__['Green2'][mesh1]:
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tComputing matrix 2 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"for depth={depth:.2e} and k={wavenumber:.2e}")
            if depth == np.infty:
                lamda_exp = np.empty(31, dtype=FLOAT_PRECISION)
                a_exp = np.empty(31, dtype=FLOAT_PRECISION)
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    mesh1.faces_centers, mesh1.faces_normals,
                    mesh2.faces_centers, mesh2.faces_areas,
                    wavenumber,         0.0,
                    self.XR, self.XZ, self.APD,
                    lamda_exp, a_exp, n_exp,
                    mesh1 is mesh2
                    )
            else:
                # Get the last computed exponential decomposition.
                a_exp, lamda_exp = next(reversed(self.exponential_decompositions.values()))
                n_exp = 31

                S2, V2 = _Green.green_2.build_matrix_2(
                    mesh1.faces_centers, mesh1.faces_normals,
                    mesh2.faces_centers, mesh2.faces_areas,
                    wavenumber, depth,
                    self.XR, self.XZ, self.APD,
                    lamda_exp, a_exp, n_exp,
                    mesh1 is mesh2
                    )

            if self.keep_matrices:
                self.__cache__['Green2'][mesh1][(mesh2, depth, wavenumber)] = (S2, V2)

        else:
            S2, V2 = self.__cache__['Green2'][mesh1][(mesh2, depth, wavenumber)]
            LOG.debug("\t" * len(_rec_depth) +
                      f"\tRetrieving stored matrix 2 of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
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
        mesh : Mesh or CollectionOfMeshes
            a mesh

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
            result.body.mesh,
            free_surface=result.free_surface,
            sea_bottom=result.sea_bottom,
            wavenumber=result.wavenumber
        )

        phi = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface, keep_details=False):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        free_surface : FreeSurface
            a meshed free surface
        keep_details : bool, optional
            if True, keep the free surface elevation in the LinearPotentialFlowResult

        Returns
        -------
        array
            the free surface elevation on each faces of the meshed free surface
        """
        fs_elevation = 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface.mesh)
        if keep_details:
            result.fs_elevation[free_surface] = fs_elevation
        return fs_elevation

