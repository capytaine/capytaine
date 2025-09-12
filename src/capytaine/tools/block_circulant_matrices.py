"""Implementation of block circulant matrices to be used for optimizing resolution with symmetries."""
# Copyright (C) 2025 Capytaine developers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
import numpy as np
from typing import List
from numpy.typing import ArrayLike
import scipy.linalg as sl

LOG = logging.getLogger(__name__)


def circular_permutation(l: List, i: int) -> List:
    return l[-i:] + l[:-i]


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
            self.nb_blocks*self.blocks[0].shape[1]
        )
        self.dtype = self.blocks[0].dtype

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        full_blocks = [np.array(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        first_row = [full_blocks[0], *(full_blocks[1:][::-1])]
        return np.block([[b for b in circular_permutation(first_row, i)]
                         for i in range(self.nb_blocks)]).astype(dtype)

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

    def block_diagonalize(self) -> "BlockDiagonalMatrix":
        if all(isinstance(b, BlockCirculantMatrix) for b in self.blocks) and self.nb_blocks == 2:
            a, b = self.blocks
            return BlockDiagonalMatrix([(a + b).block_diagonalize(), (a - b).block_diagonalize()])
        elif all(isinstance(b, np.ndarray) for b in self.blocks):
            return BlockDiagonalMatrix(np.fft.fft(np.array(self.blocks), axis=0))
        else:
            raise NotImplementedError()

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s with rhs: %s", self.__class__.__name__, self.shape, b)
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
        full_blocks = [np.array(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        return np.block([
            [full_blocks[i] if i == j else np.zeros(full_blocks[i].shape)
             for j in range(self.nb_blocks)]
            for i in range(self.nb_blocks)])

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s with rhs: %s", self.__class__.__name__, self.shape, b)
        n = self.nb_blocks
        rhs = np.split(b, n)
        res = [np.linalg.solve(Ai, bi) if isinstance(Ai, np.ndarray) else Ai.solve(bi)
               for (Ai, bi) in zip(self.blocks, rhs)]
        LOG.debug("Done")
        return np.hstack(res)


def lu_decompose(A: ArrayLike):
    if isinstance(A, np.ndarray):
        return LUDecomposedMatrix(A)
    elif isinstance(A, BlockDiagonalMatrix):
        return LUDecomposedBlockDiagonalMatrix(A)
    elif isinstance(A, BlockCirculantMatrix):
        return LUDecomposedBlockCirculantMatrix(A)
    else:
        raise NotImplementedError


class LUDecomposedMatrix:
    def __init__(self, A: np.ndarray):
        self._lu_decomp = sl.lu_factor(A)

    def solve(self, b: np.ndarray) -> np.ndarray:
        return sl.lu_solve(self._lu_decomp, b)


class LUDecomposedBlockDiagonalMatrix:
    """LU decomposition of a BlockDiagonalMatrix,
    stored as the LU decomposition of each block."""
    def __init__(self, bdm: BlockDiagonalMatrix):
        self._lu_decomp = [lu_decompose(bl) for bl in bdm.blocks]

    def solve(self, b: np.ndarray) -> np.ndarray:
        rhs = np.split(b, len(self._lu_decomp))
        res = [Ai.solve(bi) for (Ai, bi) in zip(self._lu_decomp, rhs)]
        return np.hstack(res)


class LUDecomposedBlockCirculantMatrix:
    def __init__(self, bcm: BlockCirculantMatrix):
        self._lu_decomp = lu_decompose(bcm.block_diagonalize())
        self.nb_blocks = bcm.nb_blocks

    def solve(self, b: np.ndarray) -> np.ndarray:
        n = self.nb_blocks
        b_fft = np.fft.fft(b.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        res_fft = self._lu_decomp.solve(b_fft)
        res = np.fft.ifft(res_fft.reshape((n, b.shape[0]//n)), axis=0).reshape(b.shape)
        return res
