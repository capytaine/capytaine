#!/usr/bin/env python
# coding: utf-8

import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix


################################################################################
#                       Block symmetric Toeplitz matrix                        #
################################################################################

class BlockSymmetricToeplitzMatrix(BlockMatrix):
    """A (2D) block symmetric Toeplitz matrix, stored as a list of blocks.

    Stored in the backend as a 1×N array of arrays.
    """

    # ACCESSING DATA

    def _index_grid(self):
        n = self.nb_blocks[0]
        return np.array([[abs(i-j) for i in range(n)] for j in range(n)])

    @property
    def all_blocks(self):
        return np.array([[block for block in self._stored_blocks_flat[indices]] for indices in self._index_grid()])

    @property
    def nb_blocks(self):
        return self._nb_stored_blocks[1], self._nb_stored_blocks[1]

    @property
    def block_shapes(self):
        return ([self._stored_block_shapes[0][0]]*self.nb_blocks[0],
                self._stored_block_shapes[1])

    @property
    def shape(self):
        return self._stored_block_shapes[0][0]*self.nb_blocks[0], self._stored_block_shapes[1][0]*self.nb_blocks[1]

    @property
    def block_shape(self):
        return self._stored_block_shapes[0][0], self._stored_block_shapes[1][0]  # Shape of first stored block

    @property
    def first_block_line(self):
        return self._stored_blocks_flat

    def _check_dimension(self) -> None:
        block_shape = self._stored_blocks_flat[0].shape
        for block in self._stored_blocks_flat:
            assert block.shape == block_shape  # All blocks have same shape

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        return [(i, j) for i in range(n) for j in range(n) if abs(i-j) == k]

    # DISPLAYING DATA

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        for k, block in enumerate(self._stored_blocks_flat):
            shift = np.array((global_shift[0] + k*self.block_shape[1], global_shift[1]))
            if isinstance(block, BlockMatrix):
                patches_of_block = block._patches(shift)
            elif isinstance(block, np.ndarray):
                patches_of_block = [Rectangle(shift, self.block_shape[1], self.block_shape[0], edgecolor='k', facecolor=next(self.display_color))]
            else:
                raise AttributeError()

            for i, j in self._positions_of_index(k)[1:]:
                for patch in patches_of_block:
                    local_shift = np.array(patch.get_xy()) + np.array(((j-k)*self.block_shape[1], i*self.block_shape[0]))
                    patches.append(Rectangle(local_shift, patch.get_width(), patch.get_height(), facecolor=patch.get_facecolor(), alpha=0.5))
            patches.extend(patches_of_block)

        return patches

    # TRANSFORMING DATA

    @property
    def T(self):
        transposed_blocks = np.array([[block.T for block in self._stored_blocks_flat]])
        return self.__class__(transposed_blocks)


###########################################################################
#                    Block symmetric circulant matrix                     #
###########################################################################

class BlockSymmetricCirculantMatrix(BlockSymmetricToeplitzMatrix):

    @property
    def nb_blocks(self):
        return self._nb_blocks, self._nb_blocks

    def _index_grid(self):
        line = self._baseline_grid()
        grid = [line]
        for i in range(1, self.nb_blocks[0]):
            grid.append(line[-i:] + line[:-i])
        return np.array(grid)

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        grid = self._index_grid()
        return [(i, j) for i in range(n) for j in range(n) if grid[i, j] == k]

    @property
    def first_block_line(self):
        return self._stored_blocks_flat[self._baseline_grid()]


class EvenBlockSymmetricCirculantMatrix(BlockSymmetricCirculantMatrix):
    """
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

    def __init__(self, blocks):
        super().__init__(blocks)
        self._nb_blocks = (self._nb_stored_blocks[1]-1)*2

    # ACCESSING DATA

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._stored_blocks_flat)))
        return blocks_indices[:-1] + blocks_indices[1:][::-1]


class OddBlockSymmetricCirculantMatrix(BlockSymmetricCirculantMatrix):
    """
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

    def __init__(self, blocks):
        super().__init__(blocks)
        self._nb_blocks = self._nb_stored_blocks[1]*2 - 1

    # ACCESSING DATA

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._stored_blocks_flat)))
        return blocks_indices[:] + blocks_indices[1:][::-1]

