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
        if max_rank is None and tol == 0.0:
            LOG.warning("No stopping criterion for the Adaptive Cross Approximation.")

        if max_rank is None:
            max_rank = min(nb_rows, nb_cols)//2

        A = []  # Left matrix to be assembled
        B = []  # Right matrix to be assembled
        R = np.zeros((nb_rows, nb_cols), dtype=dtype)  # Current best approximation == A @ B
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

        if tol > 0:
            LOG.warning(f"Unable to find a low rank approximation of rank lower or equal to {l+1} with tolerance {tol}.")
        return LowRankMatrix(np.array(A).T, np.array(B))

    @classmethod
    def from_rows_and_cols_functions_with_2_in_1_ACA(cls, get_row_func, get_col_func, nb_rows, nb_cols, max_rank=None, tol=0.0, dtype=np.float64):
        """Create two low rank matrices from functions running two Adaptive Cross Approximation at the same time.
        The user should provide either the `max_rank` optional argument or the `tol` optional argument to expect a useful output.

        Parameters
        ----------
        get_row_func: Function
            Function such that `get_row_func(i)` returns the `i`-th row of both full matrices.
        get_col_func: Function
            Function such that `get_col_func(j)` returns the `j`-th column of both full matrices.
        nb_rows: int
            Number of rows in both full matrices.
        nb_cols: int
            Number of columns in both full matrices.
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

        left = ([], [])  # Left matrices to be assembled
        right = ([], [])  # Right matrices to be assembled
        full0 = np.zeros((nb_rows, nb_cols), dtype=dtype)
        full1 = np.zeros((nb_rows, nb_cols), dtype=dtype)  # Current best approximations == A @ B
        available_rows = list(range(nb_rows))
        available_cols = list(range(nb_cols))

        for l in range(max_rank):

            # Peek a row
            if l == 0:
                relative_i = 0  # Or chose at random
            else:
                relative_i = int(np.argmax(np.abs(left[1][l-1][available_rows])))
            i = available_rows.pop(relative_i)
            # relative_i is the index of the row in the list of remaining rows,
            # e.g. if available_rows = [2, 7, 8] and relative_i = 2, the chosen
            # row has index 8 in the original full matrix.

            # Add the row to the approximation
            b0, b1 = get_row_func(i)
            right[0].append(b0 - full0[i, :])
            right[1].append(b1 - full1[i, :])

            # Peek a column
            relative_j = int(np.argmax(np.abs(right[1][l][available_cols])))
            j = available_cols.pop(relative_j)

            # Get the pivot
            a0, a1 = get_col_func(j)

            new_col0 = a0 - full0[:, j]
            pivot0 = new_col0[i]
            if abs(pivot0) < 1e-12:
                pivot0 = 1e-12

            new_col1 = a1 - full1[:, j]
            pivot1 = new_col1[i]
            if abs(pivot1) < 1e-12:
                pivot1 = 1e-12

            # Add the column to the approximation
            left[0].append(new_col0/pivot0)
            left[1].append(new_col1/pivot1)

            # Update the current best approximation
            increment_of_current_iteration0 = np.outer(left[0][l][:], right[0][l][:])
            full0 += increment_of_current_iteration0

            increment_of_current_iteration1 = np.outer(left[1][l][:], right[1][l][:])
            full1 += increment_of_current_iteration1

            if np.linalg.norm(increment_of_current_iteration1, 'fro') <= tol*np.linalg.norm(full1, 'fro'):
                # See Gypsilab for possible improvement of the norm computation.
                LOG.debug(f"ACA: approximation found of rank {l}")
                if l == 0:  # Edge case of the zero matrix...
                    return (LowRankMatrix(np.array(left[0]).T, np.array(right[0])),
                            LowRankMatrix(np.array(left[1]).T, np.array(right[1])))
                else:
                    return (LowRankMatrix(np.array(left[0][:-1]).T, np.array(right[0][:-1])),
                            LowRankMatrix(np.array(left[1][:-1]).T, np.array(right[1][:-1])))

        if tol > 0:
            LOG.warning(f"Unable to find a low rank approximation of rank lower or equal to {l+1} with tolerance {tol}.")
        return (LowRankMatrix(np.array(left[0]).T, np.array(right[0])),
                LowRankMatrix(np.array(left[1]).T, np.array(right[1])))

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

