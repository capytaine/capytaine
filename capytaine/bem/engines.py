#!/usr/bin/env python
# coding: utf-8
"""Definition of the methods to build influence matrices, using possibly some sparse structures."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from abc import ABC, abstractmethod

import numpy as np

from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import ReflectionSymmetricMesh, TranslationalSymmetricMesh, AxialSymmetricMesh

from capytaine.matrices import linear_solvers
from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.low_rank import LowRankMatrix, NoConvergenceOfACA
from capytaine.matrices.block_toeplitz import BlockSymmetricToeplitzMatrix, BlockToeplitzMatrix, BlockCirculantMatrix
from capytaine.tools.lru_cache import delete_first_lru_cache

LOG = logging.getLogger(__name__)


####################
#  ABSTRACT CLASS  #
####################

class MatrixEngine(ABC):
    """Abstract method to build a matrix."""

    @abstractmethod
    def build_matrices(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber, green_function):
        pass

    def build_S_matrix(self, *args, **kwargs):
        """Similar to :code:`build_matrices`, but returning only :math:`S`"""
        S, _ = self.build_matrices(*args, **kwargs)  # Could be optimized...
        return S


##################
#  BASIC ENGINE  #
##################

class BasicMatrixEngine(MatrixEngine):
    """
    Simple engine that assemble a full matrix (except for one reflection symmetry).
    Basically only calls :code:`green_function.evaluate`.

    Parameters
    ----------
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver
        (available: "direct" and "gmres", the former is the default choice)
        or by passing directly a solver function.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """

    available_linear_solvers = {'direct': linear_solvers.solve_directly,
                                'gmres': linear_solvers.solve_gmres}

    def __init__(self, *, linear_solver='direct', matrix_cache_size=1):

        if linear_solver in self.available_linear_solvers:
            self.linear_solver = self.available_linear_solvers[linear_solver]
        else:
            self.linear_solver = linear_solver

        if matrix_cache_size > 0:
            self.build_matrices = delete_first_lru_cache(maxsize=matrix_cache_size)(self.build_matrices)

        self.exportable_settings = {
            'engine': 'BasicMatrixEngine',
            'matrix_cache_size': matrix_cache_size,
            'linear_solver': str(linear_solver),
        }

    def build_matrices(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber, green_function):
        r"""Build the influence matrices between mesh1 and mesh2.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float
            position of the free surface (default: :math:`z = 0`)
        sea_bottom: float
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float
            wavenumber (default: 1.0)
        green_function: AbstractGreenFunction
            object with an "evaluate" method that computes the Green function.

        Returns
        -------
        tuple of matrix-like
            the matrices :math:`S` and :math:`K`
        """

        if (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            S_a, V_a = self.build_matrices(
                mesh1[0], mesh2[0], free_surface, sea_bottom, wavenumber,
                green_function)
            S_b, V_b = self.build_matrices(
                mesh1[0], mesh2[1], free_surface, sea_bottom, wavenumber,
                green_function)

            return BlockSymmetricToeplitzMatrix([[S_a, S_b]]), BlockSymmetricToeplitzMatrix([[V_a, V_b]])

        else:
            return green_function.evaluate(
                mesh1, mesh2, free_surface, sea_bottom, wavenumber,
            )

###################################
#  HIERARCHIAL TOEPLITZ MATRICES  #
###################################

class HierarchicalToeplitzMatrixEngine(MatrixEngine):
    """An experimental matrix engine that build a hierarchical matrix with
     some block-Toeplitz structure.

    Parameters
    ----------
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """

    def __init__(self, *, ACA_distance=8.0, ACA_tol=1e-2, matrix_cache_size=1):

        if matrix_cache_size > 0:
            self.build_matrices = delete_first_lru_cache(maxsize=matrix_cache_size)(self.build_matrices)

        self.ACA_distance = ACA_distance
        self.ACA_tol = ACA_tol

        self.linear_solver = linear_solvers.solve_gmres

        self.exportable_settings = {
            'engine': 'HierarchicalToeplitzMatrixEngine',
            'ACA_distance': ACA_distance,
            'ACA_tol': ACA_tol,
            'matrix_cache_size': matrix_cache_size,
        }

    def build_matrices(self,
                       mesh1, mesh2, *args,
                       _rec_depth=1, **kwargs):
        """Recursively builds a hierarchical matrix between mesh1 and mesh2.
        
        Same arguments as :func:`BasicMatrixEngine.build_matrices`.

        :code:`_rec_depth` keeps track of the recursion depth only for pretty log printing.
        """

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            log_entry = (
                "\t" * (_rec_depth+1) +
                "Build the S and K influence matrices between {mesh1} and {mesh2}"
                .format(mesh1=mesh1.name, mesh2=(mesh2.name if mesh2 is not mesh1 else 'itself'))
            )
        else:
            log_entry = ""  # will not be used

        green_function = args[-1]

        # Distance between the meshes (for ACA).
        distance = np.linalg.norm(mesh1.center_of_mass_of_nodes - mesh2.center_of_mass_of_nodes)

        # I) SPARSE COMPUTATION
        # I-i) BLOCK TOEPLITZ MATRIX

        if (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            LOG.debug(log_entry + " using mirror symmetry.")

            S_a, V_a = self.build_matrices(
                mesh1[0], mesh2[0], *args, **kwargs,
                _rec_depth=_rec_depth+1)
            S_b, V_b = self.build_matrices(
                mesh1[0], mesh2[1], *args, **kwargs,
                _rec_depth=_rec_depth+1)

            return BlockSymmetricToeplitzMatrix([[S_a, S_b]]), BlockSymmetricToeplitzMatrix([[V_a, V_b]])

        elif (isinstance(mesh1, TranslationalSymmetricMesh)
              and isinstance(mesh2, TranslationalSymmetricMesh)
              and np.allclose(mesh1.translation, mesh2.translation)
              and mesh1.nb_submeshes == mesh2.nb_submeshes):

            LOG.debug(log_entry + " using translational symmetry.")

            S_list, V_list = [], []
            for submesh in mesh2:
                S, V = self.build_matrices(
                    mesh1[0], submesh, *args, **kwargs,
                    _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)
            for submesh in mesh1[1:][::-1]:
                S, V = self.build_matrices(
                    submesh, mesh2[0], *args, **kwargs,
                    _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)

            return BlockToeplitzMatrix([S_list]), BlockToeplitzMatrix([V_list])

        elif (isinstance(mesh1, AxialSymmetricMesh)
              and isinstance(mesh2, AxialSymmetricMesh)
              and mesh1.axis == mesh2.axis
              and mesh1.nb_submeshes == mesh2.nb_submeshes):

            LOG.debug(log_entry + " using rotation symmetry.")

            S_line, V_line = [], []
            for submesh in mesh2[:mesh2.nb_submeshes]:
                S, V = self.build_matrices(
                    mesh1[0], submesh, *args, **kwargs,
                    _rec_depth=_rec_depth+1)
                S_line.append(S)
                V_line.append(V)

            return BlockCirculantMatrix([S_line]), BlockCirculantMatrix([V_line])

        # I-ii) LOW-RANK MATRIX WITH ACA

        elif distance > self.ACA_distance*mesh1.diameter_of_nodes or distance > self.ACA_distance*mesh2.diameter_of_nodes:

            LOG.debug(log_entry + " using ACA.")

            def get_row_func(i):
                s, v = green_function.evaluate(
                    mesh1.extract_one_face(i), mesh2,
                    *args[:-1], **kwargs
                )
                return s.flatten(), v.flatten()

            def get_col_func(j):
                s, v = green_function.evaluate(
                    mesh1, mesh2.extract_one_face(j),
                    *args[:-1], **kwargs
                )
                return s.flatten(), v.flatten()

            try:
                return LowRankMatrix.from_rows_and_cols_functions_with_multi_ACA(
                    get_row_func, get_col_func, mesh1.nb_faces, mesh2.nb_faces,
                    nb_matrices=2, id_main=1,  # Approximate V and get an approximation of S at the same time
                    tol=self.ACA_tol, dtype=np.complex128)
            except NoConvergenceOfACA:
                pass  # Continue with non sparse computation

        # II) NON-SPARSE COMPUTATIONS
        # II-i) BLOCK MATRIX

        if (isinstance(mesh1, CollectionOfMeshes)
              and isinstance(mesh2, CollectionOfMeshes)):

            LOG.debug(log_entry + " using block matrix structure.")

            S_matrix, V_matrix = [], []
            for submesh1 in mesh1:
                S_line, V_line = [], []
                for submesh2 in mesh2:
                    S, V = self.build_matrices(
                        submesh1, submesh2, *args, **kwargs,
                        _rec_depth=_rec_depth+1)

                    S_line.append(S)
                    V_line.append(V)
                S_matrix.append(S_line)
                V_matrix.append(V_line)

            return BlockMatrix(S_matrix), BlockMatrix(V_matrix)

        # II-ii) PLAIN NUMPY ARRAY

        else:
            LOG.debug(log_entry)

            S, V = green_function.evaluate(
                mesh1, mesh2, *args[:-1], **kwargs
            )
            return S, V

