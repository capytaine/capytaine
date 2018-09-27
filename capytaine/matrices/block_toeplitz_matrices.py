#!/usr/bin/env python
# coding: utf-8

import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix


class BlockSymmetricToeplitzMatrix(BlockMatrix):
    """A (2D) block symmetric Toeplitz matrix, stored as a list of blocks."""

    # ACCESSING DATA

    def _index_grid(self):
        n = self.nb_blocks[0]
        return np.array([[abs(i-j) for i in range(n)] for j in range(n)])

    @property
    def all_blocks(self):
        return np.array([[block for block in self._stored_blocks_flat[indices]] for indices in self._index_grid()])

    @property
    def block_shape(self):
        return self._stored_blocks[0][0].shape

    def _check_dimension(self) -> None:
        for block in self._stored_blocks_flat:
            assert block.shape == self.block_shape  # All blocks have same shape

    @property
    def shape(self):
        return self.nb_blocks[0]*self.block_shape[0], self.nb_blocks[1]*self.block_shape[1]

    @property
    def nb_blocks(self):
        return self._stored_blocks_flat.shape[0], self._stored_blocks_flat.shape[0]

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
        return BlockSymmetricToeplitzMatrix(transposed_blocks)


class BlockSymmetricCirculantMatrix(BlockSymmetricToeplitzMatrix):
    def __init__(self, c_blocks, size=None):
        super().__init__(c_blocks)

        n = len(self._stored_blocks_flat)
        if size is None or size == 2*n-2:
            self.is_even = True
        elif size == 2*n-1:
            self.is_even = False
        else:
            raise ValueError("Invalid size for BlockCirculantMatrix.")

    # ACCESSING DATA

    def _check_dimension(self) -> None:
        for block in self._stored_blocks_flat:
            assert block.shape == self._stored_blocks[0][0].shape  # All blocks have same shape

    @property
    def nb_blocks(self):
        if self.is_even:
            return 2*len(self._stored_blocks_flat)-2, 2*len(self._stored_blocks_flat)-2
        else:
            return 2*len(self._stored_blocks_flat)-1, 2*len(self._stored_blocks_flat)-1

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._stored_blocks_flat)))
        if self.is_even:
            base_line = blocks_indices[:-1] + blocks_indices[1:][::-1]
        else:
            base_line = blocks_indices + blocks_indices[1:][::-1]
        return base_line

    @property
    def _t_blocks(self):
        return self._stored_blocks_flat[self._baseline_grid()]

    def _index_grid(self):
        line = self._baseline_grid()
        grid = [line]
        for i in range(1, self.nb_blocks[0]):
            grid.append(line[-i:] + line[:-i])
        return np.array(grid)

    @property
    def block_shape(self):
        return self._stored_blocks_flat[0].shape

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        grid = self._index_grid()
        return [(i, j) for i in range(n) for j in range(n) if grid[i, j] == k]

    # TRANSFORMING DATA

    @property
    def T(self):
        transposed_blocks = np.array([[block.T for block in self._stored_blocks_flat]])
        return BlockSymmetricCirculantMatrix(transposed_blocks)

