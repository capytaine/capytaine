#!/usr/bin/env python
# coding: utf-8

import numpy as np
from capytaine.matrices.block_matrices import BlockMatrix


# @abc
# class SingleBlock(AbstractMatrix):
#     nb_blocks = (1, 1)

#     @property
#     def block_shapes(self):
#         return self.shape


# # class FullBlock(SingleBlock):


# class SparseZeroBlock(SingleBlock):

#     def __init__(self, shape):
#         pass


class LowRankBlock:

    ##############
    #  Creation  #
    ##############

    def __init__(self, left_matrix, right_matrix):
        self.left_matrix = left_matrix
        self.right_matrix = right_matrix
        self.shape = left_matrix.shape[0], right_matrix.shape[1]
        assert left_matrix.shape[1] == right_matrix.shape[0], "Sizes of the left and right matrices do not match."
        self.rank = left_matrix.shape[1] # == right_matrix.shape[0]

    @classmethod
    def from_full_matrix_with_SVD(cls, matrix, max_rank):
        """Create a low rank matrix from a full matrix using SVD decomposition."""
        u, s, v = np.linalg.svd(matrix)
        left_matrix = u[:, 0:max_rank] @ np.diag(s[0:max_rank])
        right_matrix = v[0:max_rank, :]
        return cls(left_matrix, right_matrix)

    @classmethod
    def from_full_matrix_with_ACA(cls, matrix, tol=1e-6, max_rank=None):
        if max_rank is None:
            max_rank = min(matrix.shape)

        A = []  # Left matrix to be assembled
        B = []  # Right matrix to be assembled
        R = np.zeros(matrix.shape)  # Current best approximation == A @ B
        available_lines = list(range(matrix.shape[0]))
        available_cols = list(range(matrix.shape[1]))

        for l in range(max_rank):

            # Peek a line
            if l == 0:
                relative_i = 0  # Or chose at random
            else:
                relative_i = np.argmax(np.abs(A[-1][available_lines]))
            i = available_lines.pop(relative_i)

            # Add the line to the approximation
            B.append(matrix[i, :] - R[i, :])

            # Peek a column
            relative_j = np.argmax(np.abs(B[-1][available_cols]))
            j = available_cols.pop(relative_j)

            # Get the pivot
            pivot = matrix[i, j] - R[i, j]
            assert pivot != 0.0

            # Take the column j for the approximation
            A.append((matrix[:, j] - R[:, j]) / pivot)

            print("i, j:", i, j)

            # Update the current best approximation
            R += np.outer(A[-1][:], B[-1][:])

            if np.linalg.norm(matrix - R, ord='fro') < tol:
                print(f"Approximation found of rank {l+1}")
                break

        return LowRankBlock(np.array(A).T, np.array(B))

    ####################
    #  Representation  #
    ####################

    def full_matrix(self):
        return self.left_matrix @ self.right_matrix

    ####################
    #  Transformation  #
    ####################

    def __matmul__(self, other):
        if isinstance(other, np.ndarray) and len(other.shape) == 1:
            return self._mul_with_vector(other)
        else:
            return NotImplemented

    def _mul_with_vector(self, other):
        return self.left_matrix @ (self.right_matrix @ other)


if __name__ == "__main__":
    # A = LowRankBlock(np.random.rand(4, 1), np.random.rand(1, 4))
    # assert A.shape == A.full_matrix().shape
    # assert np.linalg.matrix_rank(A.full_matrix()) == A.rank == 1

    n = 5
    A = np.random.rand(n, n)
    dumb_low_rank = LowRankBlock.from_full_matrix_with_SVD(A, n)
    assert np.allclose(dumb_low_rank.full_matrix() - A, 0.0)

    rank_1 = LowRankBlock.from_full_matrix_with_SVD(A, 1)
    assert np.linalg.matrix_rank(rank_1.full_matrix()) == rank_1.rank == 1

    b = np.random.rand(n)
    assert np.allclose(rank_1 @ b, rank_1.full_matrix() @ b)

    A = rank_1.full_matrix()
    print(A)
    other_rank_1 = LowRankBlock.from_full_matrix_with_ACA(A, max_rank=2)
    print(np.linalg.matrix_rank(other_rank_1.full_matrix()))
    print(A - other_rank_1.full_matrix())




