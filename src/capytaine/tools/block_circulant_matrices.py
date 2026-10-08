# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Implementation of block circulant matrices to be used for optimizing resolution with symmetries."""

import cmath
import logging
import numpy as np
from abc import ABC, abstractmethod
from functools import lru_cache, singledispatch
from typing import Any, List, Union, Sequence
from numpy.typing import NDArray
import scipy.linalg as sl

from capytaine.tools.array_backend import (
    IterableTensor, array_namespace, device, ending_dimensions_at_the_beginning, is_array,
    leading_dimensions_at_the_end, split, to_numpy,
)

LOG = logging.getLogger(__name__)


def circular_permutation(l: List, i: int) -> List:
    return l[-i:] + l[:-i]


class _BlockMatrixBackendMixin:
    """Make the array library and the device of the blocks accessible, as for an array."""

    @property
    def _namespace(self):
        return array_namespace(self.blocks[0])

    def __array_namespace__(self, *, api_version=None):
        return self._namespace

    @property
    def device(self):
        return device(self.blocks[0])


class BlockCirculantMatrix(_BlockMatrixBackendMixin):
    """Data-sparse representation of a block matrix of the following form::

        ( a  d  c  b )
        ( b  a  d  c )
        ( c  b  a  d )
        ( d  c  b  a )

    where a, b, c and d are matrices of the same shape.

    Parameters
    ----------
    blocks: Sequence of arrays, can be also an array of shape (nb_blocks, n, n, ...)
        The **first column** of blocks [a, b, c, d, ...]
        Each block should have the same shape.
    """
    def __init__(self, blocks: Union[Sequence[Any], Any]):
        self.blocks = IterableTensor(blocks)
        self.nb_blocks = len(self.blocks)
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
        full_blocks = [to_numpy(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        first_row = [full_blocks[0], *(full_blocks[1:][::-1])]
        if self.ndim >= 3:
            first_row = [leading_dimensions_at_the_end(b) for b in first_row]
            # Need to permute_dims to conform to `block` usage when the array is more than 2D
        full_matrix = np.block([[b for b in circular_permutation(first_row, i)]
                         for i in range(self.nb_blocks)])
        if self.ndim >= 3:
            full_matrix = ending_dimensions_at_the_beginning(full_matrix)
        # `self.dtype` is not used because it can be a dtype of another array library.
        return full_matrix if dtype is None else full_matrix.astype(dtype)

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
        if not (is_array(other) and other.ndim == 1):
            return NotImplemented
        xp = self._namespace
        n = other.shape[0]
        if self.nb_blocks == 2:
            a, b = self.blocks
            x1, x2 = other[:n//2], other[n//2:]
            return xp.concat([a @ x1 + b @ x2, b @ x1 + a @ x2], axis=0)
        elif self.nb_blocks == 3:
            a, b, c = self.blocks
            x1, x2, x3 = other[:n//3], other[n//3:2*n//3], other[2*n//3:]
            return xp.concat([
                a @ x1 + c @ x2 + b @ x3,
                b @ x1 + a @ x2 + c @ x3,
                c @ x1 + b @ x2 + a @ x3,
            ], axis=0)
        else:
            y = xp.zeros(other.shape, dtype=xp.result_type(self.dtype, other.dtype), device=device(other))
            blocks_indices = list(range(self.nb_blocks))
            for i, x_i in enumerate(split(other, self.nb_blocks)):
                y += xp.concat([self.blocks[j] @ x_i for j in circular_permutation(blocks_indices, i)], axis=0)
            return y

    def matvec(self, other):
        return self.__matmul__(other)

    def block_diagonalize(self) -> "BlockDiagonalMatrix":
        if self.ndim == 2 and self.nb_blocks == 2:
            a, b = self.blocks
            return BlockDiagonalMatrix([a + b, a - b])
        elif self.ndim == 2 and self.nb_blocks == 3:
            a, b, c = self.blocks
            w = cmath.exp(-2j*cmath.pi/3)
            return BlockDiagonalMatrix([
                a + b + c,
                a + w * b + w.conjugate() * c,
                a + w.conjugate() * b + w * c,
            ])
        elif self.ndim == 2 and self.nb_blocks == 4:
            a, b, c, d = self.blocks
            return BlockDiagonalMatrix([
                a + b + c + d,
                a - 1j*b - c + 1j*d,
                a - b + c - d,
                a + 1j*b - c - 1j*d,
            ])
        elif self.ndim == 2:
            xp = self._namespace
            return BlockDiagonalMatrix(xp.fft.fft(self.blocks.as_array(xp), axis=0))
        else:
            raise NotImplementedError()

    def solve(self, b):
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        xp, n = self._namespace, self.nb_blocks
        b_fft = xp.reshape(xp.fft.fft(xp.reshape(b, (n, -1)), axis=0), b.shape)
        res_fft = self.block_diagonalize().solve(b_fft)
        res = xp.reshape(xp.fft.ifft(xp.reshape(res_fft, (n, -1)), axis=0), b.shape)
        LOG.debug("Done")
        return res


class NestedBlockCirculantMatrix(_BlockMatrixBackendMixin):
    """Data-sparse representation of a block matrix of the following form::

       ( a  b | e  d | c  f )
       ( b  a | f  c | d  e )
       ( ------------------ )
       ( c  f | a  b | e  d )
       ( d  e | b  a | f  c )
       ( ------------------ )
       ( e  d | c  f | a  b )
       ( f  c | d  e | b  a )

    where a, b, c, d, e and f are matrices of the same shape,
    that is a block circulant matrix (here of size 3x3, but arbitrary outer sizes are supported),
    where diagonal blocks are 2x2 block circulant matrices and off-diagonal blocks have subblocks in common.

    Reordering the lines and columns, this matrix is equivalent to the following shape::

       ( a  e  c | b  d  f )
       ( c  a  e | f  b  d )
       ( e  c  a | d  f  b )
       ( ----------------- )
       ( b  f  d | a  c  e )
       ( d  b  f | e  a  c )
       ( f  d  b | c  e  a )

    that is a 2x2 block matrix of block circulant matrices, the block circulant matrices of the bottom row being the transpose of the block circulant matrices of the top row.

    In the 2x2x2x2 limit case, the matrix is a nested block circulant matrix::

       ( a  b | c  d )
       ( b  a | d  c )
       ( ----------- )
       ( c  d | a  b )
       ( d  c | b  a )

    The 4x4x2x2 pattern reads::

       ( a  b | g  d | e  f | c  h )
       ( b  a | h  c | f  e | d  g )
       ( -----|------|------|----- )
       ( c  h | a  b | g  d | e  f )
       ( d  g | b  a | h  c | f  e )
       ( -----|------|------|----- )
       ( e  f | c  h | a  b | g  d )
       ( f  e | d  g | b  a | h  c )
       ( -----|------|------|----- )
       ( g  d | e  f | c  h | a  b )
       ( h  c | f  e | d  g | b  a )

    Parameters
    ----------
    blocks: Sequence of arrays, can be also an array of shape (nb_blocks, n, n, ...)
        The **first column** of blocks [a, b, c, d, ...]
        Each block should have the same shape, and the number of blocks should be even.
    """
    def __init__(self, blocks: Union[Sequence[Any], Any]):
        self.blocks = IterableTensor(blocks)
        self.nb_blocks = len(self.blocks)
        assert self.nb_blocks % 2 == 0
        assert all(self.blocks[0].shape == b.shape for b in self.blocks[1:])
        assert all(self.blocks[0].dtype == b.dtype for b in self.blocks[1:])
        self.shape = (
            self.nb_blocks*self.blocks[0].shape[0],
            self.nb_blocks*self.blocks[0].shape[1],
            *self.blocks[0].shape[2:]
        )
        self.ndim = len(self.shape)
        self.dtype = self.blocks[0].dtype

    @lru_cache
    def to_BlockCirculantMatrix(self):
        """Convert to a BlockCirculantMatrix with combined macro-blocks.

        The NestedBlockCirculantMatrix with blocks ``[a, b, c, d, e, f, ...]``
        (where ``nb_blocks = 2*n``) is converted to a BlockCirculantMatrix with n macro-blocks.

        Each macro-block i (for i in 0..n-1) is a 2x2 arrangement of the original blocks:

        * For i=0: ``[[blocks[0], blocks[1]], [blocks[1], blocks[0]]]``
        * For i>0: ``[[blocks[2*i], blocks[2*n-2*i+1]], [blocks[2*i+1], blocks[2*n-2*i]]]``
        """
        n = self.nb_blocks // 2
        xp = self._namespace

        # Create the macro-blocks for the BlockCirculantMatrix
        macro_blocks = []
        for i in range(n):
            if i == 0:
                # First macro-block is a standard 2x2 circulant
                top_row = xp.concat([self.blocks[0], self.blocks[1]], axis=1)
                bottom_row = xp.concat([self.blocks[1], self.blocks[0]], axis=1)
            else:
                # Other macro-blocks follow the pattern
                idx_top_left = 2 * i
                idx_top_right = 2 * n - 2 * i + 1
                idx_bottom_left = 2 * i + 1
                idx_bottom_right = 2 * n - 2 * i

                top_row = xp.concat([self.blocks[idx_top_left], self.blocks[idx_top_right]], axis=1)
                bottom_row = xp.concat([self.blocks[idx_bottom_left], self.blocks[idx_bottom_right]], axis=1)

            macro_block = xp.concat([top_row, bottom_row], axis=0)
            macro_blocks.append(macro_block)

        # Create the outer BlockCirculantMatrix
        return BlockCirculantMatrix(macro_blocks)

    def __array__(self, dtype=None, copy=True):
        return self.to_BlockCirculantMatrix().__array__(dtype=dtype, copy=copy)

    def __matmul__(self, other):
        return self.to_BlockCirculantMatrix() @ other

    def matvec(self, other):
        return self.__matmul__(other)

    def block_diagonalize(self) -> "BlockDiagonalMatrix":
        if self.nb_blocks == 4:
            # Special case for 2x2x2x2: fully block-diagonalize both outer and inner structures
            # Given blocks [a, b, c, d], the matrix is:
            # ( a  b | c  d )
            # ( b  a | d  c )
            # ( ----------- )
            # ( c  d | a  b )
            # ( d  c | b  a )
            # First diagonalize outer 2x2: gives BlockCirculantMatrix([a+c, b+d]) and BlockCirculantMatrix([a-c, b-d])
            # Then diagonalize each inner 2x2: (a+c)±(b+d) and (a-c)±(b-d)
            a, b, c, d = self.blocks
            return BlockDiagonalMatrix([
                a + b + c + d,
                a - b + c - d,
                a + b - c - d,
                a - b - c + d
            ])
        else:
            # Drop the inner structure
            return self.to_BlockCirculantMatrix().block_diagonalize()

    def solve(self, b):
        return self.to_BlockCirculantMatrix().solve(b)


class BlockDiagonalMatrix(_BlockMatrixBackendMixin):
    """Data-sparse representation of a block matrix of the following form::

        ( a  0  0  0 )
        ( 0  b  0  0 )
        ( 0  0  c  0 )
        ( 0  0  0  d )

    where a, b, c and d are matrices of the same shape.

    Parameters
    ----------
    blocks: Sequence of arrays, can be also an array of shape (nb_blocks, n, n)
        The blocks [a, b, c, d, ...]
    """
    def __init__(self, blocks: Union[Sequence[Any], Any]):
        blocks = self.blocks = IterableTensor(blocks)
        self.nb_blocks = len(blocks)
        assert all(blocks[0].shape == b.shape for b in blocks[1:])
        self.shape = (
                sum(bl.shape[0] for bl in blocks),
                sum(bl.shape[1] for bl in blocks)
                )
        assert all(blocks[0].dtype == b.dtype for b in blocks[1:])
        self.ndim = len(self.shape)
        self.dtype = blocks[0].dtype

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        full_blocks = [to_numpy(b) for b in self.blocks]  # Transform all blocks to numpy arrays
        if self.ndim >= 3:
            full_blocks = [leading_dimensions_at_the_end(b) for b in full_blocks]
        full_matrix = np.block([
            [full_blocks[i] if i == j else np.zeros(full_blocks[i].shape, dtype=full_blocks[i].dtype)
             for j in range(self.nb_blocks)]
            for i in range(self.nb_blocks)])
        if self.ndim >= 3:
            full_matrix = ending_dimensions_at_the_beginning(full_matrix)
        return full_matrix if dtype is None else full_matrix.astype(dtype)

    def solve(self, b):
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        xp = self._namespace
        rhs = split(b, self.nb_blocks)
        # The right-hand side is always given to `solve` as a 2D array (with a single column
        # for a vector), because the behaviour of `solve` for a 1D right-hand side depends
        # on the version of the array API standard.
        res = [xp.reshape(xp.linalg.solve(Ai, xp.reshape(bi, (bi.shape[0], -1))), bi.shape)
               for (Ai, bi) in zip(self.blocks, rhs)]
        LOG.debug("Done")
        return xp.concat(res, axis=0)


MatrixLike = Union[np.ndarray, BlockDiagonalMatrix, NestedBlockCirculantMatrix, BlockCirculantMatrix]


class AbstractLUDecomposedMatrix(ABC):
    """Base class of the LU decompositions of matrices.

    New matrix types (e.g. matrices stored with another array library) can be
    supported by the linear solvers of Capytaine by registering a function
    returning a subclass of this class with :func:`lu_decompose`.
    """
    shape: tuple
    dtype: object

    @property
    def _namespace(self):
        return self.__array_namespace__()

    @abstractmethod
    def __array_namespace__(self, *, api_version=None):
        """Array library of the matrix that has been decomposed, as for an array of the array API standard."""

    @property
    @abstractmethod
    def device(self):
        """Device of the matrix that has been decomposed, as for an array of the array API standard."""

    @abstractmethod
    def solve(self, b):
        """Solve the linear system with the decomposed matrix as left-hand side and `b` as right-hand side."""


@singledispatch
def lu_decompose(A, *, overwrite_a: bool = False) -> AbstractLUDecomposedMatrix:
    """Compute the LU decomposition of `A`.

    This is a :func:`functools.singledispatch` function: the implementation for
    a new type of matrix can be added with ``lu_decompose.register(MyMatrix, my_function)``,
    where ``my_function(A, *, overwrite_a=False)`` returns an :class:`AbstractLUDecomposedMatrix`.
    """
    raise NotImplementedError(f"No LU decomposition registered for {type(A)}")


class LUDecomposedMatrix(AbstractLUDecomposedMatrix):
    def __init__(self, A: NDArray, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  A.__class__.__name__, A.shape)
        self._lu_decomp = sl.lu_factor(A, overwrite_a=overwrite_a)
        self.shape = A.shape
        self.dtype = A.dtype
        self._array_namespace = array_namespace(A)
        self._device = device(A)

    def __array_namespace__(self, *, api_version=None):
        return self._array_namespace

    @property
    def device(self):
        return self._device

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        return sl.lu_solve(self._lu_decomp, b)


class LUDecomposedBlockDiagonalMatrix(AbstractLUDecomposedMatrix):
    """LU decomposition of a BlockDiagonalMatrix,
    stored as the LU decomposition of each block."""
    def __init__(self, bdm: BlockDiagonalMatrix, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  bdm.__class__.__name__, bdm.shape)
        self._lu_decomp = [lu_decompose(bl, overwrite_a=overwrite_a) for bl in bdm.blocks]
        self.shape = bdm.shape
        self.nb_blocks = bdm.nb_blocks
        self.dtype = bdm.dtype

    def __array_namespace__(self, *, api_version=None):
        return self._lu_decomp[0].__array_namespace__()

    @property
    def device(self):
        return self._lu_decomp[0].device

    def solve(self, b):
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        xp = self._namespace
        rhs = split(b, self.nb_blocks)
        res = [Ai.solve(bi) for (Ai, bi) in zip(self._lu_decomp, rhs)]
        return xp.concat(res, axis=0)


class LUDecomposedBlockCirculantMatrix(AbstractLUDecomposedMatrix):
    def __init__(self, bcm: BlockCirculantMatrix, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  bcm.__class__.__name__, bcm.shape)
        self._lu_decomp = lu_decompose(bcm.block_diagonalize(), overwrite_a=overwrite_a)
        self.shape = bcm.shape
        self.nb_blocks = bcm.nb_blocks
        self.dtype = bcm.dtype

    def __array_namespace__(self, *, api_version=None):
        return self._lu_decomp.__array_namespace__()

    @property
    def device(self):
        return self._lu_decomp.device

    def solve(self, b):
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        xp, n = self._namespace, self.nb_blocks
        b_fft = xp.reshape(xp.fft.fft(xp.reshape(b, (n, -1)), axis=0), b.shape)
        res_fft = self._lu_decomp.solve(b_fft)
        return xp.reshape(xp.fft.ifft(xp.reshape(res_fft, (n, -1)), axis=0), b.shape)


@lu_decompose.register(np.ndarray)
def _lu_decompose_ndarray(A, *, overwrite_a: bool = False):
    return LUDecomposedMatrix(A, overwrite_a=overwrite_a)


@lu_decompose.register(BlockDiagonalMatrix)
def _lu_decompose_block_diagonal(A, *, overwrite_a: bool = False):
    return LUDecomposedBlockDiagonalMatrix(A, overwrite_a=overwrite_a)


@lu_decompose.register(BlockCirculantMatrix)
def _lu_decompose_block_circulant(A, *, overwrite_a: bool = False):
    return LUDecomposedBlockCirculantMatrix(A, overwrite_a=overwrite_a)


@lu_decompose.register(NestedBlockCirculantMatrix)
def _lu_decompose_nested_block_circulant(A, *, overwrite_a: bool = False):
    return LUDecomposedBlockCirculantMatrix(A.to_BlockCirculantMatrix(), overwrite_a=overwrite_a)


def has_been_lu_decomposed(A):
    return isinstance(A, AbstractLUDecomposedMatrix)
