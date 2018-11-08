#!/usr/bin/env python
# coding: utf-8

import logging
from typing import Tuple, List, Iterable
from functools import lru_cache

import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix

LOG = logging.getLogger(__name__)


################################################################################
#                       Block symmetric Toeplitz matrix                        #
################################################################################

class BlockSymmetricToeplitzMatrix(BlockMatrix):
    """A (2D) block symmetric Toeplitz matrix, stored as a list of blocks.
    All blocks should have the same shape.

    Stored in the backend as a 1×N array of arrays.
    """

    def _compute_shape(self):
        # The full shape is found by multiplying the shape of the blocks. All of them have the same shape.
        return self._stored_block_shapes[0][0]*self.nb_blocks[0], self._stored_block_shapes[1][0]*self.nb_blocks[1]

    def _check_dimensions_of_blocks(self) -> bool:
        for block in self._stored_blocks[0, 1:]:
            if not block.shape == self.block_shape:  # All blocks have same shape
                return False
        return True

    # ACCESSING DATA

    def _index_grid(self) -> np.ndarray:
        """Helper function to find the positions at which the blocks appears in the full matrix.

        Example of output::

            [[1, 2, 3, 4, 5],
             [2, 1, 2, 3, 4],
             [3, 2, 1, 2, 3],
             [4, 3, 2, 1, 2],
             [5, 4, 3, 2, 1]]
        """
        n = self.nb_blocks[0]
        return np.array([[abs(i-j) for i in range(n)] for j in range(n)])

    @property
    def first_block_line(self):
        """The blocks on the first line of blocks in the matrix."""
        return self._stored_blocks[0, :]

    @property
    def all_blocks(self) -> np.ndarray:
        """The matrix of matrices as if the block Toeplitz structure was not used."""
        return np.array([[block for block in self.first_block_line[indices]] for indices in self._index_grid()])

    @property
    def nb_blocks(self) -> Tuple[int, int]:
        """The number of blocks in each direction."""
        return self._stored_nb_blocks[1], self._stored_nb_blocks[1]

    @property
    def block_shapes(self):
        """The shapes of the blocks composing the block matrix.
        Actually, they should be all the same."""
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                self._stored_block_shapes[1])

    @property
    def block_shape(self):
        """The shape of any block."""
        return self._stored_block_shapes[0][0], self._stored_block_shapes[1][0]

    def _block_indices_of(self, k: int) -> List[Tuple[int, int]]:
        """The block indices at which the block k from the first line can also be found.
        TODO: Optimize.
        """
        n = self.nb_blocks[0]
        return [(i, j) for i in range(n) for j in range(n) if abs(i-j) == k]

    def _positions_of(self, k: int, global_frame=(0, 0)) -> List[Tuple[int, int]]:
        """The positions in the full matrix at which the block k from the first line can also be found."""
        shape = self.block_shape
        return [(global_frame[0] + i*shape[0], global_frame[1] + j*shape[1]) for i, j in self._block_indices_of(k)]

    def _stored_block_positions(self, global_frame=(0, 0)) -> Iterable[List[Tuple[int, int]]]:
        """The position of each blocks in the matrix.

        Example::

            AABB
            AABB  ->  list(matrix._stored_block_positions) = [[(0,0), (2, 2)], [(0, 2), (2, 0)]]
            BBAA
            BBAA
        """
        return (self._positions_of(k, global_frame=global_frame) for k in range(len(self.first_block_line)))

    # TRANSFORMING DATA

    @lru_cache(maxsize=16)
    def _circulant_super_matrix(self):
        return EvenBlockSymmetricCirculantMatrix(self._stored_blocks,
                                                 _stored_block_shapes=self._stored_block_shapes,
                                                 check=False)

    def matvec(self, other):
        A = self._circulant_super_matrix()
        b = np.concatenate([other, np.zeros(A.shape[1] - self.shape[1])])
        return (A @ b)[:self.shape[0]]

    def rmatvec(self, other):
        A = self._circulant_super_matrix()
        b = np.concatenate([other, np.zeros(A.shape[0] - self.shape[0])])
        return (A.rmatvec(b))[:self.shape[1]]

    @property
    def T(self):
        """Transpose the matrix."""
        return self._apply_unary_op(lambda x: x.T)


###########################################################################
#                    Block symmetric circulant matrix                     #
###########################################################################

class _AbstractBlockSymmetricCirculantMatrix(BlockSymmetricToeplitzMatrix):
    """Should not be instantiated. Just here to factor some common code between the two classes below."""

    # ACCESSING DATA

    def _index_grid(self):
        line = self._baseline_grid()
        grid = [line]
        for i in range(1, self.nb_blocks[0]):
            grid.append(line[-i:] + line[:-i])
        return np.array(grid)

    def _block_indices_of(self, k):
        n = self.nb_blocks[0]
        grid = self._index_grid()
        return [(i, j) for i in range(n) for j in range(n) if grid[i, j] == k]

    @property
    def first_block_line(self):
        return self._stored_blocks[0, self._baseline_grid()]

    @property
    def nb_blocks(self):
        return self._nb_blocks, self._nb_blocks

    @property
    def block_shapes(self):
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                [self._stored_block_shapes[1][0]]*self.nb_blocks[1])

    # TRANSFORMING DATA

    @lru_cache(maxsize=16)
    def block_diagonalize(self):
        stacked_blocks = np.array([block.full_matrix() if not isinstance(block, np.ndarray) else block
                       for block in self.first_block_line])
        blocks_of_diagonalization = np.fft.fft(stacked_blocks, axis=0)
        return blocks_of_diagonalization

    def matvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        fft_of_vector = np.fft.fft(np.reshape(other, (self.nb_blocks[0], self.block_shape[1], 1)), axis=0)
        fft_of_result = self.block_diagonalize() @ fft_of_vector
        result = np.fft.ifft(fft_of_result, axis=0).reshape(self.shape[0])
        return np.asarray(result, dtype=other.dtype)

    def rmatvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        other = np.conjugate(other)
        fft_of_vector = np.fft.fft(np.reshape(other, (self.nb_blocks[0], 1, self.block_shape[0])), axis=0)
        fft_of_result = fft_of_vector @ self.block_diagonalize()
        result = np.fft.ifft(fft_of_result, axis=0).reshape(self.shape[1])
        return np.asarray(result, dtype=other.dtype)


class EvenBlockSymmetricCirculantMatrix(_AbstractBlockSymmetricCirculantMatrix):
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
    """

    def __init__(self, blocks, **kwargs):
        blocks = np.asarray(blocks)
        self._nb_blocks = (blocks.shape[1] - 1) * 2
        super().__init__(blocks, **kwargs)

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._stored_blocks[0, :])))
        return blocks_indices[:-1] + blocks_indices[1:][::-1]


class OddBlockSymmetricCirculantMatrix(_AbstractBlockSymmetricCirculantMatrix):
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
    """

    def __init__(self, blocks, **kwargs):
        blocks = np.asarray(blocks)
        self._nb_blocks = (blocks.shape[1]) * 2 - 1
        super().__init__(blocks, **kwargs)

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._stored_blocks[0, :])))
        return blocks_indices[:] + blocks_indices[1:][::-1]

