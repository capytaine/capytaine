#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np

LOG = logging.getLogger(__name__)


class LowRankMatrix:

    ##############
    #  Creation  #
    ##############

    def __init__(self, left_matrix, right_matrix):
        self.left_matrix = left_matrix
        self.right_matrix = right_matrix
        self.shape = left_matrix.shape[0], right_matrix.shape[1]
        assert left_matrix.shape[1] == right_matrix.shape[0], "Sizes of the left and right matrices do not match."
        self.rank = left_matrix.shape[1]  # == right_matrix.shape[0]
        assert left_matrix.dtype == right_matrix.dtype
        self.dtype = left_matrix.dtype

    @classmethod
    def from_full_matrix_with_SVD(cls, full_matrix, max_rank):
        """Create a low rank matrix from a full matrix using Singular Value Decomposition."""
        u, s, v = np.linalg.svd(full_matrix)
        left_matrix = u[:, 0:max_rank] @ np.diag(s[0:max_rank])
        right_matrix = v[0:max_rank, :]
        return cls(left_matrix, right_matrix)

    @classmethod
    def from_full_matrix_with_ACA(cls, full_matrix, max_rank=None, tol=1e-6):
        """Create a low rank matrix from a full matrix using Adaptive Cross Approximation"""
        def get_row(i):
            return full_matrix[i, :]

        def get_col(j):
            return full_matrix[:, j]

        return cls.from_rows_and_cols_functions_with_ACA(
            get_row, get_col, full_matrix.shape[0], full_matrix.shape[1], max_rank=max_rank, tol=tol
        )

    @classmethod
    def from_function_with_ACA(cls, func, nb_rows, nb_cols, max_rank=None, tol=1e-6):
        def get_row(i):
            return np.asarray([func(i, j) for j in range(nb_cols)])

        def get_col(j):
            return np.asarray([func(i, j) for i in range(nb_rows)])

        return cls.from_rows_and_cols_functions_with_ACA(
            get_row, get_col, nb_rows, nb_cols, max_rank=max_rank, tol=tol
        )

    @classmethod
    def from_rows_and_cols_functions_with_ACA(cls, get_row_func, get_col_func, nb_rows, nb_cols, max_rank=None, tol=1e-6):
        """Create a low rank matrix from functions using Adaptive Cross Approximation"""
        if max_rank is None:
            max_rank = min(nb_rows, nb_cols)//2

        A = []  # Left matrix to be assembled
        B = []  # Right matrix to be assembled
        R = np.zeros((nb_rows, nb_cols))  # Current best approximation == A @ B
        available_rows = list(range(nb_rows))
        available_cols = list(range(nb_cols))

        for l in range(max_rank):

            # Peek a row
            if l == 0:
                relative_i = 0  # Or chose at random
            else:
                relative_i = np.argmax(np.abs(A[l-1][available_rows]))
            i = available_rows.pop(relative_i)
            # relative_i is the index of the row in the list of remaining rows,
            # e.g. if available_rows = [2, 7, 8] and relative_i = 2, the chosen
            # row has index 8 in the original full matrix.

            # Add the row to the approximation
            B.append(get_row_func(i) - R[i, :])

            # Peek a column
            relative_j = np.argmax(np.abs(B[l][available_cols]))
            j = available_cols.pop(relative_j)

            # Get the pivot
            new_col = get_col_func(j) - R[:, j]
            pivot = new_col[i]
            assert pivot != 0.0

            # Add the column to the approximation
            A.append(new_col/pivot)

            # Update the current best approximation
            increment_of_current_iteration = np.outer(A[l][:], B[l][:])
            R += increment_of_current_iteration

            if np.linalg.norm(increment_of_current_iteration, 'fro')/np.linalg.norm(R, 'fro') < tol:
                # See Gypsilab for possible improvement of the norm computation.
                LOG.debug(f"ACA: approximation found of rank {l}")
                return LowRankMatrix(np.array(A[:-1]).T, np.array(B[:-1]))  # Drop the last iteration

        LOG.warning(f"Unable to find a suitable low rank approximation of rank lower or equal to {l+1}.")
        return LowRankMatrix(np.array(A).T, np.array(B))

    ####################
    #  Representation  #
    ####################

    def full_matrix(self):
        return self.left_matrix @ self.right_matrix

    ####################
    #  Transformation  #
    ####################

    def recompress(self, tol=None, new_rank=None):
        """From Gipsylab, hmxQRSVD.m"""
        if new_rank is None:
            new_rank = self.rank
        QA, RA = np.linalg.qr(self.left_matrix)
        QB, RB = np.linalg.qr(self.right_matrix.T)
        U, S, V = np.linalg.svd(RA @ RB.T)
        if tol is not None:
            new_rank = np.count_nonzero(S/S[0] >= tol)
        A = QA @ (U[:, :new_rank] @ np.diag(S[:new_rank]))
        B = QB @ V[:, :new_rank]
        return LowRankMatrix(A, B.T)

    def __add__(self, other):
        if isinstance(other, LowRankMatrix):
            new_left = np.concatenate([self.left_matrix, other.left_matrix], axis=1)
            new_right = np.concatenate([self.right_matrix, other.right_matrix], axis=0)
            return LowRankMatrix(new_left, new_right).recompress(new_rank=min(self.rank, other.rank))
        else:
            return NotImplemented

    def __matmul__(self, other):
        if isinstance(other, np.ndarray) and len(other.shape) == 1:
            return self._mul_with_vector(other)
        else:
            return NotImplemented

    def _mul_with_vector(self, other):
        return self.left_matrix @ (self.right_matrix @ other)

