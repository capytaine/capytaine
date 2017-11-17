#!/usr/bin/env python
# coding: utf-8

import logging
from itertools import product

import numpy as np

LOG = logging.getLogger(__name__)


class BlockToeplitzMatrix:
    """A symmetric block Toeplitz matrix stored as a list of matrices."""

    def __init__(self, blocks):
        """
        Parameters
        ----------
        blocks: list of square matrices
            the blocks of the first row (or the first column) of the block matrix.
            they should be square matrices of the same size and the same type.
        """

        self.blocks = []
        self.dtype = blocks[0].dtype

        for block in blocks:
            if isinstance(block, BlockToeplitzMatrix):
                # Recursive block matrices not implemented yet.
                self.blocks.append(block.full_matrix())
            else:
                self.blocks.append(block)

        for block in self.blocks:
            assert len(block.shape) == 2
            assert block.shape[0] == self.block_size
            assert block.shape[1] == self.block_size
            assert block.dtype == self.dtype

    @property
    def shape(self):
        return self.nb_blocks*self.block_size, self.nb_blocks*self.block_size

    @property
    def nb_blocks(self):
        return len(self.blocks)

    @property
    def block_size(self):
        return self.blocks[0].shape[0]

    def __eq__(self, other):
        if isinstance(other, BlockToeplitzMatrix):
            # Return true if the content is the same
            # even if the decomposition is not the same
            return self.full_matrix() == other.full_matrix()
        elif isinstance(other, np.ndarray):
            return self.full_matrix() == other
        else:
            raise NotImplemented

    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            new_blocks = []
            for i in range(self.nb_blocks):
                new_blocks.append(self.blocks[i] + other)
            return BlockToeplitzMatrix(new_blocks)
        elif isinstance(other, BlockToeplitzMatrix):
            # Keep the symmetric block Toeplitz structure
            new_blocks = []
            for i in range(self.nb_blocks):
                new_blocks.append(self.blocks[i] + other.blocks[i])
            return BlockToeplitzMatrix(new_blocks)
        else:
            # Lose the symmetric block Toeplitz structure
            return self.full_matrix() + other

    def __radd__(self, other):
        # Addition is commutative
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __neg__(self):
        new_blocks = []
        for i in range(self.nb_blocks):
            new_blocks.append(- self.blocks[i])
        return BlockToeplitzMatrix(new_blocks)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            new_blocks = []
            for i in range(self.nb_blocks):
                new_blocks.append(self.blocks[i] * other)
            return BlockToeplitzMatrix(new_blocks)
        elif isinstance(other, np.ndarray):
            return self.full_matrix() * other
        else:
            raise NotImplemented

    def __rmul__(self, other):
        # Multiplication is commutative
        return self.__mul__(other)

    def __matmul__(self, other):
        if isinstance(other, np.ndarray):
            if self.nb_blocks*self.block_size != other.shape[0]:
                raise Exception("Size of the matrices does not match!")
            result = np.zeros(other.shape, dtype=self.dtype)
            for i, j in product(range(self.nb_blocks), repeat=2):
                result[i*self.block_size:(i+1)*self.block_size] += \
                    self.blocks[abs(i-j)] @ other[j*self.block_size:(j+1)*self.block_size]
            return result
        else:
            raise NotImplemented

    def __rmatmul__(self, other):
        if isinstance(other, np.ndarray):
            return other @ self.full_matrix()
        else:
            raise NotImplemented

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            new_blocks = []
            for i in range(self.nb_blocks):
                new_blocks.append(self.blocks[i] / other)
            return BlockToeplitzMatrix(new_blocks)
        else:
            raise NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            new_blocks = []
            for i in range(self.nb_blocks):
                new_blocks.append(other / self.blocks[i])
            return BlockToeplitzMatrix(new_blocks)
        else:
            raise NotImplemented

    def full_matrix(self):
        """Return the matrix as an usual array not using the symmetry."""
        full_matrix = np.empty(self.shape, dtype=self.dtype)
        for i in range(self.nb_blocks):
            for j in range(self.nb_blocks):
                full_matrix[i*self.block_size:(i+1)*self.block_size,
                            j*self.block_size:(j+1)*self.block_size] = self.blocks[abs(j-i)]
        return full_matrix


def block_Toeplitz_identity(nb_blocks, block_size, **kwargs):
    """Return the identity matrix as a block Toeplitz matrix of specified size."""
    return BlockToeplitzMatrix(
        [np.identity(block_size, **kwargs)] +
        [np.zeros((block_size, block_size), **kwargs) for _ in range(nb_blocks - 1)]
    )


class BlockCirculantMatrix(BlockToeplitzMatrix):
    """A symmetric block circulant matrix stored as a list of matrices.
    """

    def __init__(self, blocks):
        """
        Parameters
        ----------
        blocks: list of square matrices
            half of the blocks of the first row (or the first column) of the block matrix.
            they should be square matrices of the same size and the same type.
        """
        BlockToeplitzMatrix.__init__(self, blocks + blocks[-2:0:-1])

    # @property
    # def nb_blocks(self):
    #     return 2*(len(self.blocks)-1)

    # def full_matrix(self):
    #     """Return the matrix as an usual array not using the symmetry."""
    #     full_matrix = np.empty(self.shape, dtype=self.dtype)
    #     for i in range(self.nb_blocks):
    #         for j in range(self.nb_blocks):
    #             if abs(i-j) < self.nb_blocks//2 + 1:
    #                 i_block = abs(i-j)
    #             else:
    #                 i_block = self.nb_blocks - abs(i-j)
    #             full_matrix[i*self.block_size:(i+1)*self.block_size,
    #                         j*self.block_size:(j+1)*self.block_size] = self.blocks[i_block]
    #     return full_matrix


# def block_circulant_identity(nb_blocks, block_size, **kwargs):
#     """Return the identity matrix as a block Circulant matrix of specified size."""
#     return BlockCirculantMatrix(
#         [np.identity(block_size, **kwargs)] +
#         [np.zeros((block_size, block_size), **kwargs) for _ in range(nb_blocks//2 - 1)]
#     )


def solve(A, b):
    """Solve the linear system Ax = b"""
    if isinstance(A, BlockCirculantMatrix):
        LOG.debug("\tSolve linear system %ix%i BlockCirculantMatrix (block size: %i)",
                  A.nb_blocks, A.nb_blocks, A.block_size)
        AA = np.stack(A.blocks)
        AAt = np.fft.fft(AA, axis=0)
        b = np.reshape(b, (A.nb_blocks, A.block_size))
        bt = np.fft.fft(b, axis=0)
        xt = solve(AAt, bt)
        x = np.fft.ifft(xt, axis=0)
        return x.reshape(A.nb_blocks*A.block_size)

    elif isinstance(A, BlockToeplitzMatrix):
        if A.nb_blocks == 2:
            LOG.debug("\tSolve system of 2x2 BlockToeplitzMatrix (block size: %i)", A.block_size)
            A1, A2 = A.blocks
            b1, b2 = b[:len(b)//2], b[len(b)//2:]
            x_plus = solve(A1 + A2, b1 + b2)
            x_minus = solve(A1 - A2, b1 - b2)
            return np.concatenate([x_plus + x_minus, x_plus - x_minus])/2

        else:
            LOG.debug("\tSolve linear system %ix%i BlockToeplitzMatrix (block size: %i)", A.nb_blocks, A.nb_blocks, A.block_size)
            # Not implemented yet
            return np.linalg.solve(A.full_matrix(), b)

    elif isinstance(A, np.ndarray):
        LOG.debug("\tSolve linear system (size: %i) with numpy.", A.shape[0])
        return np.linalg.solve(A, b)

    else:
        raise ValueError(f"Unrecognized type of {A} in solve")
