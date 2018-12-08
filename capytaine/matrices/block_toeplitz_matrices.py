#!/usr/bin/env python
# coding: utf-8
"""Block symmetric Toeplitz matrices and block symmetric circulant matrices
to be used in hierarchical matrices.
"""

import logging
from typing import Tuple, List, Iterable

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

    # INITIALIZATION

    def _compute_shape(self) -> Tuple[int, int]:
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
    def first_block_line(self) -> np.ndarray:
        """1D array of the the blocks of the first line."""
        return self._stored_blocks[0, :]

    @property
    def all_blocks(self) -> np.ndarray:
        """The matrix of blocks as if the block Toeplitz structure was not used."""
        all_blocks = np.empty(self.nb_blocks, dtype=np.object)
        all_blocks[:, :] = [[block for block in self.first_block_line[indices]] for indices in self._index_grid()]
        return all_blocks

    @property
    def nb_blocks(self) -> Tuple[int, int]:
        """The number of blocks in each direction."""
        return self._stored_nb_blocks[1], self._stored_nb_blocks[1]

    @property
    def block_shapes(self) -> Tuple[List[int], List[int]]:
        """The shapes of the blocks composing the block matrix.
        Actually, they should be all the same."""
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                self._stored_block_shapes[1])

    @property
    def block_shape(self) -> Tuple[int, int]:
        """The shape of any of the blocks."""
        return self._stored_block_shapes[0][0], self._stored_block_shapes[1][0]

    def _block_indices_of(self, k: int) -> List[Tuple[int, int]]:
        """The block indices at which the block k from the first line can also be found."""
        n = self.nb_blocks[0]
        # TODO: Optimize?
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

    def matvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Product of {self} with vector of shape {other.shape}")
        if not hasattr(self, 'circulant_super_matrix'):
            self.circulant_super_matrix = EvenBlockSymmetricCirculantMatrix(
                self._stored_blocks,
                _stored_block_shapes=self._stored_block_shapes,
                check=False)
        A = self.circulant_super_matrix
        b = np.concatenate([other, np.zeros(A.shape[1] - self.shape[1])])
        return (A @ b)[:self.shape[0]]

    def rmatvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Product of vector of shape {other.shape} with {self}")
        if other.ndim == 2 and other.shape[0] == 1:  # Actually a 1×N matrix
            other = other[0, :]
        if not hasattr(self, 'circulant_super_matrix'):
            self.circulant_super_matrix = EvenBlockSymmetricCirculantMatrix(
                self._stored_blocks,
                _stored_block_shapes=self._stored_block_shapes,
                check=False)
        A = self.circulant_super_matrix
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

    def _index_grid(self) -> np.ndarray:
        line = self._baseline_grid()
        grid = [line]
        for i in range(1, self.nb_blocks[0]):
            grid.append(line[-i:] + line[:-i])
        return np.array(grid)

    def _block_indices_of(self, k) -> List[Tuple[int, int]]:
        n = self.nb_blocks[0]
        grid = self._index_grid()
        return [(i, j) for i in range(n) for j in range(n) if grid[i, j] == k]

    @property
    def first_block_line(self) -> np.ndarray:
        return self._stored_blocks[0, self._baseline_grid()]

    @property
    def nb_blocks(self) -> Tuple[int, int]:
        return self._nb_blocks, self._nb_blocks

    @property
    def block_shapes(self) -> Tuple[List[int], List[int]]:
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                [self._stored_block_shapes[1][0]]*self.nb_blocks[1])

    # TRANSFORMING DATA

    def block_diagonalize(self):
        """Returns an array of matrices"""
        if not hasattr(self, 'block_diagonalization'):
            if all(isinstance(matrix, BlockMatrix) for matrix in self._stored_blocks[0, :]):
                self.block_diagonalization = BlockMatrix.fft_of_list(*self.first_block_line)
            else:
                stacked_blocks = np.empty((self.nb_blocks[0],) + self.block_shape, dtype=self.dtype)
                for i, block in enumerate(self.first_block_line):
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
        except TypeError:  # Or do the same thing with list comprehension.
            fft_of_result = np.array([block @ vec for block, vec in zip(blocks_of_diagonalization, fft_of_vector)])
        result = np.fft.ifft(fft_of_result, axis=0).reshape(self.shape[0])
        if self.dtype == np.complexfloating or other.dtype == np.complexfloating:
            return np.asarray(result)
        else:
            return np.asarray(np.real(result))

    def rmatvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        other = np.conjugate(other)
        fft_of_vector = np.fft.fft(np.reshape(other, (self.nb_blocks[0], 1, self.block_shape[0])), axis=0)
        blocks_of_diagonalization = self.block_diagonalize()
        try:  # Try to run it as vectorized numpy arrays.
            fft_of_result = fft_of_vector @ blocks_of_diagonalization
        except TypeError:  # Or do the same thing with list comprehension.
            fft_of_result = np.array([block.rmatvec(vec) for block, vec in zip(blocks_of_diagonalization, fft_of_vector)])
        result = np.fft.ifft(fft_of_result, axis=0).reshape(self.shape[1])
        if self.dtype == np.complexfloating or other.dtype == np.complexfloating:
            return np.asarray(result)
        else:
            return np.asarray(np.real(result))


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
        self._nb_blocks = (len(blocks[0]) - 1) * 2
        super().__init__(blocks, **kwargs)

    def _baseline_grid(self) -> List[int]:
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
        self._nb_blocks = len(blocks[0]) * 2 - 1
        super().__init__(blocks, **kwargs)

    def _baseline_grid(self) -> List[int]:
        blocks_indices = list(range(len(self._stored_blocks[0, :])))
        return blocks_indices[:] + blocks_indices[1:][::-1]

