#!/usr/bin/env python
# coding: utf-8
"""This modules implements block Toeplitz matrices to be used in Hierarchical Toeplitz matrices.

The module also contains several special cases such as block symmetric Toeplitz matrices and block circulant matrices.
All classes inherits from the BlockMatrix class.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from typing import Tuple, List, Set, Iterable

import numpy as np

from capytaine.matrices.block import BlockMatrix

LOG = logging.getLogger(__name__)


################################################################################
#                            Block Toeplitz matrix                             #
################################################################################

class BlockToeplitzMatrix(BlockMatrix):
    """A (2D) block Toeplitz matrix, stored as a list of blocks.
    All blocks should have the same shape.

    Stored in the backend as a 1×(2N-1) array of arrays."""

    # INITIALIZATION

    def _compute_shape(self) -> Tuple[int, int]:
        # The full shape is found by multiplying the shape of the blocks. All of them have the same shape.
        return (self._stored_block_shapes[0][0]*self.nb_blocks[0],
                self._stored_block_shapes[1][0]*self.nb_blocks[1])

    def _compute_nb_blocks(self) -> Tuple[int, int]:
        """Will be overridden by subclasses."""
        assert self._stored_nb_blocks[1] % 2 == 1, "Expecting an odd number of blocks to build a Toeplitz matrix"
        n = (self._stored_nb_blocks[1]+1)//2
        return n, n

    def _check_dimensions_of_blocks(self) -> bool:
        for block in self._stored_blocks[0, :]:
            if not block.shape == self.block_shape:  # All blocks have same shape
                return False
        return True

    # ACCESSING DATA

    @property
    def block_shapes(self) -> Tuple[List[int], List[int]]:
        """The shapes of the blocks composing the block matrix.
        Actually, they should be all the same."""
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                [self._stored_block_shapes[1][0]]*self.nb_blocks[1])

    @property
    def block_shape(self) -> Tuple[int, int]:
        """The shape of any of the blocks."""
        return self._stored_block_shapes[0][0], self._stored_block_shapes[1][0]

    def _block_indices_of(self, k: int) -> Set[Tuple[int, int]]:
        """The block indices at which the stored block k can be found in the full matrix of size n.
        Will be overridden by subclasses."""
        n = self.nb_blocks[0]

        if k < n:
            i, j = 0, k  # Upper triangle
        elif n <= k < 2*n:
            i, j = 2*n-1-k, 0  # Lower triangle
        else:
            raise AttributeError

        indices = set()
        while i < n and j < n:
            indices.add((i, j))
            i, j = i+1, j+1  # Going along the diagonal

        return indices

    @property
    def all_blocks(self):
        all_blocks = np.empty(self.nb_blocks, dtype=object)
        for k in range(self._stored_nb_blocks[1]):
            for i, j in self._block_indices_of(k):
                all_blocks[i, j] = self._stored_blocks[0, k]
        return all_blocks

    def _positions_of(self, k: int, global_frame=(0, 0)) -> List[Tuple[int, int]]:
        """The positions in the full matrix at which the block k from the first line can also be found."""
        shape = self.block_shape
        return sorted([(global_frame[0] + i*shape[0], global_frame[1] + j*shape[1])
                       for i, j in self._block_indices_of(k)])

    def _stored_block_positions(self, global_frame=(0, 0)) -> Iterable[List[Tuple[int, int]]]:
        """The position of each blocks in the matrix.

        Example::

            AABB
            AABB  ->  list(matrix._stored_block_positions) = [[(0,0), (2, 2)], [(0, 2), (2, 0)]]
            BBAA
            BBAA
        """
        return (self._positions_of(k, global_frame=global_frame) for k in range(self._stored_nb_blocks[1]))

    # # TRANSFORMING DATA

    @property
    def circulant_super_matrix(self):
        if not hasattr(self, '_circulant_super_matrix'):
            self._circulant_super_matrix = BlockCirculantMatrix(
                self._stored_blocks,
                _stored_block_shapes=self._stored_block_shapes,
                check=False)
        return self._circulant_super_matrix

    def matvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Product of {self} with vector of shape {other.shape}")
        A = self.circulant_super_matrix
        b = np.concatenate([other, np.zeros(A.shape[1] - self.shape[1])])
        return (A @ b)[:self.shape[0]]

    def rmatvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Product of vector of shape {other.shape} with {self}")
        if other.ndim == 2 and other.shape[0] == 1:  # Actually a 1×N matrix
            other = other[0, :]
        A = self.circulant_super_matrix
        b = np.concatenate([other, np.zeros(A.shape[0] - self.shape[0])])
        return (A.rmatvec(b))[:self.shape[1]]


################################################################################
#                       Block symmetric Toeplitz matrix                        #
################################################################################

class BlockSymmetricToeplitzMatrix(BlockToeplitzMatrix):
    """A (2D) block symmetric Toeplitz matrix, stored as a list of blocks.
    All blocks should have the same shape.

    Stored in the backend as a 1×N array of arrays."""

    def _compute_nb_blocks(self) -> Tuple[int, int]:
        n = self._stored_nb_blocks[1]
        return n, n

    def _block_indices_of(self, k: int) -> Set[Tuple[int, int]]:
        """The block indices at which the stored block k can be found in the full matrix."""
        n = self.nb_blocks[0]
        assert k < n

        if k == 0:
            return BlockToeplitzMatrix._block_indices_of(self, 0)
        else:
            return BlockToeplitzMatrix._block_indices_of(self, k).union(
                BlockToeplitzMatrix._block_indices_of(self, 2*n-k-1))

    @property
    def circulant_super_matrix(self):
        if not hasattr(self, '_circulant_super_matrix'):
            self._circulant_super_matrix = EvenBlockSymmetricCirculantMatrix(
                self._stored_blocks,
                _stored_block_shapes=self._stored_block_shapes,
                check=False)
        return self._circulant_super_matrix


################################################################################
#                            Block circulant matrix                            #
################################################################################

class BlockCirculantMatrix(BlockToeplitzMatrix):
    """A (2D) block circulant matrix, stored as a list of blocks.
    All blocks should have the same shape.

    Stored in the backend as a 1×N array of arrays."""

    def _compute_nb_blocks(self) -> Tuple[int, int]:
        n = self._stored_nb_blocks[1]
        return n, n

    def _block_indices_of(self, k: int) -> Set[Tuple[int, int]]:
        """The block indices at which the stored block k can be found in the full matrix."""
        n = self.nb_blocks[0]
        assert k < n

        if k == 0:
            return BlockToeplitzMatrix._block_indices_of(self, 0)
        else:
            return BlockToeplitzMatrix._block_indices_of(self, k).union(
                BlockToeplitzMatrix._block_indices_of(self, n+k-1))

    # LINEAR SYSTEMS

    def block_diagonalize(self):
        """Returns an array of matrices"""
        if not hasattr(self, 'block_diagonalization'):
            if all(isinstance(matrix, BlockMatrix) for matrix in self._stored_blocks[0, :]):
                self.block_diagonalization = BlockMatrix.fft_of_list(*self.all_blocks[:, 0])
            else:
                stacked_blocks = np.empty((self.nb_blocks[1],) + self.block_shape, dtype=self.dtype)
                for i, block in enumerate(self.all_blocks[:, 0]):
                    stacked_blocks[i] = block.full_matrix() if not isinstance(block, np.ndarray) else block
                self.block_diagonalization =  np.fft.fft(stacked_blocks, axis=0)
        return self.block_diagonalization

    def matvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Product of {self} with vector of shape {other.shape}")
        fft_of_vector = np.fft.fft(np.reshape(other, (self.nb_blocks[0], self.block_shape[1], 1)), axis=0)
        blocks_of_diagonalization = self.block_diagonalize()
        try:  # Try to run it as vectorized numpy arrays.
            fft_of_result = blocks_of_diagonalization @ fft_of_vector
        # When the above fails, numpy 1.15 returns a TypeError, whereas numpy 1.16 returns a ValueError.
        except (TypeError, ValueError):  # Or do the same thing with list comprehension.
            fft_of_result = np.array([block @ vec for block, vec in zip(blocks_of_diagonalization, fft_of_vector)])
        result = np.fft.ifft(fft_of_result, axis=0).reshape(self.shape[0])
        if self.dtype == complex or other.dtype == complex:
            return np.asarray(result)
        else:
            return np.asarray(np.real(result))

    def rmatvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        other = np.conjugate(other)
        fft_of_vector = np.fft.ifft(np.reshape(other, (self.nb_blocks[0], 1, self.block_shape[0])), axis=0)
        blocks_of_diagonalization = self.block_diagonalize()
        try:  # Try to run it as vectorized numpy arrays.
            fft_of_result = fft_of_vector @ blocks_of_diagonalization
        # When the above fails, numpy 1.15 returns a TypeError, whereas numpy 1.16 returns a ValueError.
        except (TypeError, ValueError):
            # Instead we do the same thing with list comprehension.
            fft_of_result = np.array(
                [block.rmatvec(vec.flatten()) for block, vec in zip(blocks_of_diagonalization, fft_of_vector)]
            )
        result = np.fft.fft(fft_of_result, axis=0).reshape(self.shape[1])
        if self.dtype == complex or other.dtype == complex:
            return np.asarray(result)
        else:
            return np.asarray(np.real(result))


###########################################################################
#                    Block symmetric circulant matrix                     #
###########################################################################

class EvenBlockSymmetricCirculantMatrix(BlockCirculantMatrix, BlockSymmetricToeplitzMatrix):
    """A block symmetric circulant matrix, with an even number of blocks.

    Examples::

        ABCB
        BABC
        CBAB
        BCBA

        ABCDCB
        BABCDB
        CBABCD
        DCBABC
        CDCBAB
        BCDCBA

    Stored in the backend as a 1×(N/2+1) array of arrays."""

    def _compute_nb_blocks(self) -> Tuple[int, int]:
        """The number of blocks in the full matrix."""
        n = (self._stored_nb_blocks[1] - 1)*2
        return n, n

    def _block_indices_of(self, k: int) -> List[Tuple[int, int]]:
        n = self.nb_blocks[0]
        assert k < n/2 + 1
        if k == 0:
            return BlockToeplitzMatrix._block_indices_of(self, 0)
        else:
            return (BlockToeplitzMatrix._block_indices_of(self, k) |
                BlockToeplitzMatrix._block_indices_of(self, n+k-1) |
                BlockToeplitzMatrix._block_indices_of(self, n-k)   |
                BlockToeplitzMatrix._block_indices_of(self, 2*n-k-1)
                    )


class OddBlockSymmetricCirculantMatrix(BlockCirculantMatrix, BlockSymmetricToeplitzMatrix):
    """A block symmetric circulant matrix, with an odd number of blocks.

    Examples::

        ABCCB
        BABCC
        CBABC
        CCBAB
        BCCBA

        ABCDDCB
        BABCDDB
        CBABCDD
        DCBABCD
        DDCBABC
        CDDCBAB
        BCDDCBA

    Stored in the backend as a 1×(N+1)/2 array of arrays."""

    def _compute_nb_blocks(self) -> Tuple[int, int]:
        n = self._stored_nb_blocks[1]*2 - 1
        return n, n

    def _block_indices_of(self, k: int) -> List[Tuple[int, int]]:
        n = self.nb_blocks[0]
        assert k < (n+1)/2
        if k == 0:
            return BlockToeplitzMatrix._block_indices_of(self, 0)
        else:
            return (BlockToeplitzMatrix._block_indices_of(self, k) |
                BlockToeplitzMatrix._block_indices_of(self, n+k-1) |
                BlockToeplitzMatrix._block_indices_of(self, n-k)   |
                BlockToeplitzMatrix._block_indices_of(self, 2*n-k-1)
                    )
