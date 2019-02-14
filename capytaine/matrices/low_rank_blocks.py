#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np

LOG = logging.getLogger(__name__)


class LowRankMatrix:
    """Matrix defined as the tensor product of two small matrices.

    Parameters
    ----------
    left_matrix: np.ndarray
    right_matrix: np.ndarray

    Attributes
    ----------
    shape: pair of ints
        The shape of the full matrix.
    rank: int
        The rank of the full matrix.
    """

    ndim = 2

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
    def from_full_matrix_with_ACA(cls, full_matrix, max_rank=None, tol=0.0):
        """Create a low rank matrix from a full matrix using Adaptive Cross Approximation.
        The user should provide either the `max_rank` optional argument or the `tol` optional argument to expect a useful output.

        Parameters
        ----------
        full_matrix: numpy array
            The matrix that will be approximated.
        nb_rows: int
        nb_cols: int
        max_rank: int, optional
        tol: float, optional
        dtype: numpy.dtype, optional
            See `from_rows_and_cols_functions_with_ACA`

        Returns
        -------
        LowRankMatrix
        """
        def get_row(i):
            return full_matrix[i, :]

        def get_col(j):
            return full_matrix[:, j]

        return cls.from_rows_and_cols_functions_with_ACA(
            get_row, get_col, full_matrix.shape[0], full_matrix.shape[1], max_rank=max_rank, tol=tol, dtype=full_matrix.dtype
        )

    @classmethod
    def from_function_with_ACA(cls, func, nb_rows, nb_cols, max_rank=None, tol=0.0, dtype=np.float64):
        """Create a low rank matrix from a function using Adaptive Cross Approximation.
        The user should provide either the `max_rank` optional argument or the `tol` optional argument to expect a useful output.

        Parameters
        ----------
        func: Function
            Function such that `func(i, j)` returns the value of the (i, j) entry of the full matrix.
        nb_rows: int
        nb_cols: int
        max_rank: int, optional
        tol: float, optional
        dtype: numpy.dtype, optional
            See `from_rows_and_cols_functions_with_ACA`

        Returns
        -------
        LowRankMatrix
        """
        def get_row(i):
            return np.asarray([func(i, j) for j in range(nb_cols)])

        def get_col(j):
            return np.asarray([func(i, j) for i in range(nb_rows)])

        return cls.from_rows_and_cols_functions_with_ACA(
            get_row, get_col, nb_rows, nb_cols, max_rank=max_rank, tol=tol, dtype=dtype
        )

    @classmethod
    def from_rows_and_cols_functions_with_ACA(cls, get_row_func, get_col_func, nb_rows, nb_cols, max_rank=None, tol=0.0, dtype=np.float64):
        """Create a low rank matrix from functions using Adaptive Cross Approximation.
        The user should provide either the `max_rank` optional argument or the `tol` optional argument to expect a useful output.

        Parameters
        ----------
        get_row_func: Function
            Function such that `get_row_func(i)` returns the `i`-th row of the full matrix.
        get_col_func: Function
            Function such that `get_col_func(j)` returns the `j`-th column of the full matrix.
        nb_rows: int
            Number of rows in the full matrix.
        nb_cols: int
            Number of columns in the full matrix.
        max_rank: int, optional
            The maximum rank allowed for the output low rank matrix.
            The default value is half the size of the full matrix, that is no gain in storage space.
        tol: float, optional
            The tolerance on the relative error (default: 0).
            If the Frobenius norm of the increment is lower than the tolerance, the iteration stops.
            If the tolerance is set to 0, the resulting matrix will have the maximum rank defined by `max_rank`.
        dtype: numpy.dtype, optional
            The type of data in the low rank matrix (default: float64).

        Returns
        -------
        LowRankMatrix
        """
        # Actually use the multi matrix ACA below.
        # Need some wrapping and unwrapping...
        def get_row(i):
            return [get_row_func(i)]

        def get_col(j):
            return [get_col_func(j)]

        return cls.from_rows_and_cols_functions_multi_ACA(
            get_row, get_col, nb_rows, nb_cols,
            nb_matrices=1, id_main=0,
            max_rank=max_rank, tol=tol, dtype=dtype
        )[0]

    @classmethod
    def from_rows_and_cols_functions_multi_ACA(cls, get_row_func, get_col_func, nb_rows, nb_cols,
                                               nb_matrices=1, id_main=0,
                                               max_rank=None, tol=0.0, dtype=np.float64):
        """Create two low rank matrices from functions running two Adaptive Cross Approximation at the same time.
        The user should provide either the `max_rank` optional argument or the `tol` optional argument to expect a useful output.

        Parameters
        ----------
        get_row_func: Function
            Function such that `get_row_func(i)` returns the `i`-th row of all full matrices.
        get_col_func: Function
            Function such that `get_col_func(j)` returns the `j`-th column of all full matrices.
        nb_rows: int
            Number of rows in all full matrices.
        nb_cols: int
            Number of columns in all full matrices.
        nb_matrices: int, optional
            The number of matrices approximated at the same time
        id_main: int, optional
            The matrix used primarily in the ACA
        max_rank: int, optional
            The maximum rank allowed for both output low rank matrices.
            The default value is half the size of the full matrix, that is no gain in storage space.
        tol: float, optional
            The tolerance on the relative error (default: 0).
            If the Frobenius norm of the increment is lower than the tolerance, the iteration stops.
            If the tolerance is set to 0, the resulting matrix will have the maximum rank defined by `max_rank`.
        dtype: numpy.dtype, optional
            The type of data in both low rank matrices (default: float64).

        Returns
        -------
        Tuple[LowRankMatrix, LowRankMatrix]
        """
        if max_rank is None and tol == 0.0:
            LOG.warning("No stopping criterion for the Adaptive Cross Approximation.")

        if max_rank is None:
            max_rank = min(nb_rows, nb_cols)//2

        # Initialize working matrices
        left = np.zeros((nb_matrices, nb_rows, max_rank), dtype=dtype)
        right = np.zeros((nb_matrices, max_rank, nb_cols), dtype=dtype)

        squared_norm_of_full_matrix = 0.0

        # List of indices of unused entries
        available_rows = list(range(nb_rows))
        available_cols = list(range(nb_cols))

        for l in range(max_rank):
            # Peek a row
            if l == 0:
                relative_i = 0  # Or chose at random
            else:
                relative_i = int(np.argmax(np.abs(left[id_main, available_rows, l - 1])))
            i = available_rows.pop(relative_i)
            # relative_i is the index of the row in the list of remaining rows,
            # e.g. if available_rows = [2, 7, 8] and relative_i = 2, the chosen
            # row has index 8 in the original full matrix.

            # Add the row to the approximation
            one_row = get_row_func(i)
            for id_mat in range(nb_matrices):
                right[id_mat, l, :] = one_row[id_mat] - left[id_mat, i, :l] @ right[id_mat, :l, :]

            # Peek a column
            relative_j = int(np.argmax(np.abs(right[id_main, l, available_cols])))
            j = available_cols.pop(relative_j)
            one_col = get_col_func(j)

            # Add the column to the approximation in all matrices
            for id_mat in range(nb_matrices):
                new_col = one_col[id_mat] - left[id_mat, :, :l] @ right[id_mat, :l, j]
                pivot = new_col[i]
                if abs(pivot) < 1e-12:
                    pivot = 1e-12
                left[id_mat, :, l] = new_col/pivot

            # Update norm of the full matrix
            squared_norm_of_increment = (
                    (left[id_main, :, l] @ left[id_main, :, l]) *
                    (right[id_main, l, :] @ right[id_main, l, :])
            )

            crossed_terms = (
                    (np.conj(left[id_main, :, l].T) @ left[id_main, :, :l]) @
                    (np.conj(right[id_main, l, :]) @ right[id_main, :l, :].T)
            )
            squared_norm_of_full_matrix += squared_norm_of_increment + 2*np.real(crossed_terms)

            if squared_norm_of_increment <= tol**2*squared_norm_of_full_matrix:
                LOG.debug(f"ACA: approximation found of rank {l}")

                if l == 0:  # Edge case of the zero matrix...
                    l = 1  # ... actually returns a rank 1 matrix.

                return [LowRankMatrix(left[id_mat, :, :l], right[id_mat, :l, :]) for id_mat in range(nb_matrices)]

        if tol > 0:
            LOG.warning(f"Unable to find a low rank approximation"
                        f"of rank lower or equal to {max_rank} with tolerance {tol}.")
        return [LowRankMatrix(left[id_mat, :, :], right[id_mat, :, :]) for id_mat in range(nb_matrices)]

    ####################
    #  Representation  #
    ####################

    def full_matrix(self):
        return self.left_matrix @ self.right_matrix

    @property
    def stored_data_size(self):
        return np.product(self.left_matrix.shape) + np.product(self.right_matrix.shape)

    @property
    def sparcity(self):
        return self.stored_data_size/np.product(self.shape)

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

    def __neg__(self):
        return LowRankMatrix(-self.left_matrix, self.right_matrix)

    def __sub__(self, other):
        return self + (-other)

    def __truediv__(self, other):
        from numbers import Number
        if isinstance(other, Number):
            return LowRankMatrix(self.left_matrix, self.right_matrix/other)
        else:
            return NotImplemented

    def __matmul__(self, other):
        if isinstance(other, np.ndarray) and len(other.shape) == 1:
            return self._mul_with_vector(other)
        else:
            return NotImplemented

    def _mul_with_vector(self, other):
        return self.left_matrix @ (self.right_matrix @ other)

    def astype(self, dtype):
        return LowRankMatrix(self.left_matrix.astype(dtype), self.right_matrix.astype(dtype))

