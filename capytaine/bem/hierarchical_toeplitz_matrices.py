#!/usr/bin/env python
# coding: utf-8
"""This module implements a decorator that recursively assemble a pair of Hierarchical Toeplitz matrices."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from functools import wraps

import numpy as np

from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import ReflectionSymmetricMesh, TranslationalSymmetricMesh, AxialSymmetricMesh

from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.low_rank import LowRankMatrix
from capytaine.matrices.block_toeplitz import BlockSymmetricToeplitzMatrix, BlockToeplitzMatrix, BlockCirculantMatrix


LOG = logging.getLogger(__name__)


def hierarchical_toeplitz_matrices(build_matrices, ACA_distance=8.0, ACA_tol=1e-2, dtype=np.float64):
    """Decorator for the matrix building functions.

    Parameters
    ----------
    build_matrices: function
        Function that takes as argument two meshes and several other parameters
        and returns a pair of influence matrices.
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    dtype: numpy.dtype, optional
        The type of the data in the matrix (typically float64 or complex128).

    Returns
    -------
    function
        A similar function that returns two hierarchical Toeplitz matrices.
    """

    if logging.getLogger().isEnabledFor(logging.DEBUG):
        # Read the docstring of the function and use the first line for logging.
        function_description_for_logging = build_matrices.__doc__.splitlines()[0] \
            .replace("mesh1", "{mesh1}").replace("mesh2", "{mesh2}")
    else:
        function_description_for_logging = ""  # irrelevant

    @wraps(build_matrices)  # Is this decorator really necessary?
    def build_hierarchical_toeplitz_matrix(mesh1, mesh2, *args, _rec_depth=1, **kwargs):
        """Assemble hierarchical Toeplitz matrices.

        The method is basically an ugly multiple dispatch on the kind of mesh.
        For hierarchical structures, the method is called recursively on all of the sub-bodies.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        *args, **kwargs
            Other arguments, passed to the actual evaluation of the coefficients
        _rec_depth: int, optional
            internal parameter: recursion accumulator, used only for pretty logging

        Returns
        -------
        pair of matrix-like objects
            influence matrices
        """

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            log_entry = "\t" * (_rec_depth+1) + function_description_for_logging.format(
                mesh1=mesh1.name, mesh2=(mesh2.name if mesh2 is not mesh1 else 'itself')
            )
        else:
            log_entry = ""  # irrelevant

        # Distance between the meshes (for ACA).
        distance = np.linalg.norm(mesh1.center_of_mass_of_nodes - mesh2.center_of_mass_of_nodes)

        # I) SPARSE COMPUTATION

        if (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            LOG.debug(log_entry + " using mirror symmetry.")

            S_a, V_a = build_hierarchical_toeplitz_matrix(
                mesh1[0], mesh2[0], *args, **kwargs,
                _rec_depth=_rec_depth+1)
            S_b, V_b = build_hierarchical_toeplitz_matrix(
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
                S, V = build_hierarchical_toeplitz_matrix(
                    mesh1[0], submesh, *args, **kwargs,
                    _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)
            for submesh in mesh1[1:][::-1]:
                S, V = build_hierarchical_toeplitz_matrix(
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
                S, V = build_hierarchical_toeplitz_matrix(
                    mesh1[0], submesh, *args, **kwargs,
                    _rec_depth=_rec_depth+1)
                S_line.append(S)
                V_line.append(V)

            return BlockCirculantMatrix([S_line]), BlockCirculantMatrix([V_line])

        elif distance > ACA_distance*mesh1.diameter_of_nodes or distance > ACA_distance*mesh2.diameter_of_nodes:
            # Low-rank matrix computed with Adaptive Cross Approximation.

            LOG.debug(log_entry + " using ACA.")

            def get_row_func(i):
                s, v = build_matrices(mesh1.extract_one_face(i), mesh2, *args, **kwargs)
                return s.flatten(), v.flatten()

            def get_col_func(j):
                s, v = build_matrices(mesh1, mesh2.extract_one_face(j), *args, **kwargs)
                return s.flatten(), v.flatten()

            return LowRankMatrix.from_rows_and_cols_functions_with_multi_ACA(
                get_row_func, get_col_func, mesh1.nb_faces, mesh2.nb_faces,
                nb_matrices=2, id_main=1,  # Approximate V and get an approximation of S at the same time
                tol=ACA_tol, dtype=dtype)

        # II) NON-SPARSE COMPUTATIONS

        elif (isinstance(mesh1, CollectionOfMeshes)
              and isinstance(mesh2, CollectionOfMeshes)):
            # Recursively build a block matrix

            LOG.debug(log_entry + " using block matrix structure.")

            S_matrix, V_matrix = [], []
            for submesh1 in mesh1:
                S_line, V_line = [], []
                for submesh2 in mesh2:
                    S, V = build_hierarchical_toeplitz_matrix(
                        submesh1, submesh2, *args, **kwargs,
                        _rec_depth=_rec_depth+1)

                    S_line.append(S)
                    V_line.append(V)
                S_matrix.append(S_line)
                V_matrix.append(V_line)

            return BlockMatrix(S_matrix), BlockMatrix(V_matrix)

        else:
            # Actual evaluation of coefficients using the Green function.

            LOG.debug(log_entry)

            return build_matrices(mesh1, mesh2, *args, **kwargs)

    return build_hierarchical_toeplitz_matrix
