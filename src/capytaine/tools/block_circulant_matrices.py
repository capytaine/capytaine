"""Implementation of block circulant matrices to be used for optimizing resolution with symmetries."""
# Copyright (C) 2025 Capytaine developers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
import numpy as np
from typing import List, Union
from numpy.typing import NDArray, ArrayLike
import scipy.linalg as sl

LOG = logging.getLogger(__name__)


def circular_permutation(l: List, i: int) -> List:
    return l[-i:] + l[:-i]


def leading_dimensions_at_the_end(a):
    """Transform an array of shape (n, m, ...) into (..., n, m).
    Invert of `leading_dimensions_at_the_end`"""
    return np.moveaxis(a, [0, 1], [-2, -1])


def ending_dimensions_at_the_beginning(a):
    """Transform an array of shape (..., n, m) into (n, m, ...).
    Invert of `leading_dimensions_at_the_end`"""
    return np.moveaxis(a, [-2, -1], [0, 1])


class BlockCirculantMatrix:
    """Data-sparse representation of a block matrix of the following form

        ( a  d  c  b )
        ( b  a  d  c )
        ( c  b  a  d )
        ( d  c  b  a )

    where a, b, c and d are matrices of the same shape.

    Parameters
    ----------
    blocks: iterable of matrix-like
        The **first column** of blocks [a, b, c, d, ...]
        Each block should have the same shape.
    """
    def __init__(self, blocks: List[ArrayLike]):
        self.blocks = [b for b in blocks]
        self.nb_blocks = len(blocks)
        assert all(self.blocks[0].shape == b.shape for b in self.blocks[1:])
        assert all(self.blocks[0].dtype == b.dtype for b in self.blocks[1:])
        self.shape = (
            self.nb_blocks*self.blocks[0].shape[0],
            self.nb_blocks*self.blocks[0].shape[1],
            *self.blocks[0].shape[2:]
        )
        self.ndim = len(self.shape)
        self.dtype = self.blocks[0].dtype

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        if dtype is None:
            dtype = self.dtype
        full_blocks = [np.array(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        first_row = [full_blocks[0], *(full_blocks[1:][::-1])]
        if self.ndim >= 3:
            first_row = [leading_dimensions_at_the_end(b) for b in first_row]
            # Need to permute_dims to conform to `block` usage when the array is more than 2D
        full_matrix = np.block([[b for b in circular_permutation(first_row, i)]
                         for i in range(self.nb_blocks)]).astype(dtype)
        if self.ndim >= 3:
            full_matrix = ending_dimensions_at_the_beginning(full_matrix)
        return full_matrix

    def __add__(self, other):
        if isinstance(other, BlockCirculantMatrix) and self.shape == other.shape:
            return BlockCirculantMatrix([a + b for (a, b) in zip(self.blocks, other.blocks)])
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, BlockCirculantMatrix) and self.shape == other.shape:
            return BlockCirculantMatrix([a - b for (a, b) in zip(self.blocks, other.blocks)])
        else:
            return NotImplemented

    def __matmul__(self, other):
        if self.nb_blocks == 2 and isinstance(other, np.ndarray) and other.ndim == 1:
            a, b = self.blocks
            x1, x2 = other[:len(other)//2], other[len(other)//2:]
            return np.concatenate([a @ x1 + b @ x2, b @ x1 + a @ x2], axis=0)
        else:
            return NotImplemented

    def block_diagonalize(self) -> "BlockDiagonalMatrix":
        if self.ndim == 2 and all(isinstance(b, BlockCirculantMatrix) for b in self.blocks) and self.nb_blocks == 2:
            a, b = self.blocks
            return BlockDiagonalMatrix([a + b, a - b])
        elif self.ndim == 2 and all(isinstance(b, np.ndarray) for b in self.blocks):
            return BlockDiagonalMatrix(np.fft.fft(np.array(self.blocks), axis=0))
        else:
            raise NotImplementedError()

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        n = self.nb_blocks
        b_fft = np.fft.fft(b.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        res_fft = self.block_diagonalize().solve(b_fft)
        res = np.fft.ifft(res_fft.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        LOG.debug("Done")
        return res


class BlockDiagonalMatrix:
    """Data-sparse representation of a block matrix of the following form

        ( a  0  0  0 )
        ( 0  b  0  0 )
        ( 0  0  c  0 )
        ( 0  0  0  d )

    where a, b, c and d are matrices of the same shape.

    Parameters
    ----------
    blocks: iterable of matrix-like
        The blocks [a, b, c, d, ...]
    """
    def __init__(self, blocks):
        self.blocks = [b for b in blocks]
        self.nb_blocks = len(blocks)
        assert all(blocks[0].shape == b.shape for b in blocks[1:])
        self.shape = (
                sum(bl.shape[0] for bl in blocks),
                sum(bl.shape[1] for bl in blocks)
                )
        assert all(blocks[0].dtype == b.dtype for b in blocks[1:])
        self.dtype = blocks[0].dtype

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        if dtype is None:
            dtype = self.dtype
        full_blocks = [np.array(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        if self.ndim >= 3:
            full_blocks = [leading_dimensions_at_the_end(b) for b in full_blocks]
        full_matrix = np.block([
            [full_blocks[i] if i == j else np.zeros(full_blocks[i].shape)
             for j in range(self.nb_blocks)]
            for i in range(self.nb_blocks)])
        if self.ndim >= 3:
            full_matrix = ending_dimensions_at_the_beginning(full_matrix)
        return full_matrix

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        n = self.nb_blocks
        rhs = np.split(b, n)
        res = [np.linalg.solve(Ai, bi) if isinstance(Ai, np.ndarray) else Ai.solve(bi)
               for (Ai, bi) in zip(self.blocks, rhs)]
        LOG.debug("Done")
        return np.hstack(res)


MatrixLike = Union[np.ndarray, BlockDiagonalMatrix, BlockCirculantMatrix]


def lu_decompose(A: MatrixLike, *, overwrite_a : bool = False):
    if isinstance(A, np.ndarray):
        return LUDecomposedMatrix(A, overwrite_a=overwrite_a)
    elif isinstance(A, BlockDiagonalMatrix):
        return LUDecomposedBlockDiagonalMatrix(A, overwrite_a=overwrite_a)
    elif isinstance(A, BlockCirculantMatrix):
        return LUDecomposedBlockCirculantMatrix(A, overwrite_a=overwrite_a)
    else:
        raise NotImplementedError()


class LUDecomposedMatrix:
    def __init__(self, A: NDArray, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  A.__class__.__name__, A.shape)
        self._lu_decomp = sl.lu_factor(A, overwrite_a=overwrite_a)
        self.shape = A.shape
        self.dtype = A.dtype

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        return sl.lu_solve(self._lu_decomp, b)


class LUDecomposedBlockDiagonalMatrix:
    """LU decomposition of a BlockDiagonalMatrix,
    stored as the LU decomposition of each block."""
    def __init__(self, bdm: BlockDiagonalMatrix, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  bdm.__class__.__name__, bdm.shape)
        self._lu_decomp = [lu_decompose(bl, overwrite_a=overwrite_a) for bl in bdm.blocks]
        self.shape = bdm.shape
        self.nb_blocks = bdm.nb_blocks
        self.dtype = bdm.dtype

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        rhs = np.split(b, self.nb_blocks)
        res = [Ai.solve(bi) for (Ai, bi) in zip(self._lu_decomp, rhs)]
        return np.hstack(res)


class LUDecomposedBlockCirculantMatrix:
    def __init__(self, bcm: BlockCirculantMatrix, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  bcm.__class__.__name__, bcm.shape)
        self._lu_decomp = lu_decompose(bcm.block_diagonalize(), overwrite_a=overwrite_a)
        self.shape = bcm.shape
        self.nb_blocks = bcm.nb_blocks
        self.dtype = bcm.dtype

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        n = self.nb_blocks
        b_fft = np.fft.fft(b.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        res_fft = self._lu_decomp.solve(b_fft)
        res = np.fft.ifft(res_fft.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        return res


LUDecomposedMatrixLike = Union[LUDecomposedMatrix, LUDecomposedBlockDiagonalMatrix, LUDecomposedBlockCirculantMatrix]


def has_been_lu_decomposed(A):
    # Python 3.8 does not support isinstance(A, LUDecomposedMatrixLike)
    return isinstance(A, (LUDecomposedMatrix, LUDecomposedBlockDiagonalMatrix, LUDecomposedBlockCirculantMatrix))
