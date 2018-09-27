#!/usr/bin/env python
# coding: utf-8

import numpy as np
from capytaine.matrices.block_matrices import BlockMatrix


class BlockSymmetricToeplitzMatrix(BlockMatrix):
    def __init__(self, t_blocks):
        # TODO: Check dimensionality
        self._t_blocks = np.asanyarray(t_blocks)

    @property
    def nb_blocks(self):
        return self._t_blocks.shape[0], self._t_blocks.shape[0]

    @property
    def block_shape(self):
        return self._t_blocks[0].shape

    @property
    def shape(self):
        return self.nb_blocks[0]*self.block_shape[0], self.nb_blocks[1]*self.block_shape[1]

    @property
    def T(self):
        transposed_blocks = [block.T for block in self._t_blocks]
        return BlockSymmetricToeplitzMatrix(transposed_blocks)

    @property
    def blocks(self):
        return [[block for block in self._t_blocks[indices]] for indices in self._index_grid()]

    def _index_grid(self):
        n = self.nb_blocks[0]
        return np.array([[abs(i-j) for i in range(n)] for j in range(n)])

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        return [(i, j) for i in range(n) for j in range(n) if abs(i-j) == k]

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        for k, block in enumerate(self._t_blocks):
            shift = np.array((global_shift[0] + k*self.block_shape[1], global_shift[1]))
            if isinstance(block, BlockMatrix):
                patches_of_block = block._patches(shift)
            elif isinstance(block, np.ndarray):
                patches_of_block = [Rectangle(shift, self.block_shape[1], self.block_shape[0], edgecolor='k', facecolor=next(self.COLORS))]
            else:
                raise AttributeError()

            for i, j in self._positions_of_index(k)[1:]:
                for patch in patches_of_block:
                    local_shift = np.array(patch.get_xy()) + np.array(((j-k)*self.block_shape[1], i*self.block_shape[0]))
                    patches.append(Rectangle(local_shift, patch.get_width(), patch.get_height(), facecolor=patch.get_facecolor(), alpha=0.5))
            patches.extend(patches_of_block)

        return patches


class BlockSymmetricCirculantMatrix(BlockSymmetricToeplitzMatrix):
    def __init__(self, c_blocks, size=None):
        # TODO: Check dimensionality
        self._c_blocks = np.asanyarray(c_blocks)
        if size is None or size == 2*len(c_blocks)-1:
            self.is_even = True
        elif size == 2*len(c_blocks)-2:
            self.is_even = False
        else:
            raise ValueError("Invalid size for BlockCirculantMatrix.")

    @property
    def nb_blocks(self):
        if self.is_even:
            return 2*len(self._c_blocks)-2, 2*len(self._c_blocks)-2
        else:
            return 2*len(self._c_blocks)-1, 2*len(self._c_blocks)-1

    def _baseline_grid(self):
        blocks_indices = list(range(len(self._c_blocks)))
        if self.is_even:
            base_line = blocks_indices[:-1] + blocks_indices[1:][::-1]
        else:
            base_line = blocks_indices + blocks_indices[1:][::-1]
        return base_line

    @property
    def _t_blocks(self):
        return self._c_blocks[self._baseline_grid()]

    def _index_grid(self):
        line = self._baseline_grid()
        grid = [line]
        for i in range(1, self.nb_blocks[0]):
            grid.append(line[-i:] + line[:-i])
        return np.array(grid)

    @property
    def block_shape(self):
        return self._c_blocks[0].shape

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        grid = self._index_grid()
        return [(i, j) for i in range(n) for j in range(n) if grid[i, j] == k]

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        for k, block in enumerate(self._c_blocks):
            shift = np.array((global_shift[0] + k*self.block_shape[1], global_shift[1]))
            if isinstance(block, BlockMatrix):
                patches_of_block = block._patches(shift)
            elif isinstance(block, np.ndarray):
                patches_of_block = [Rectangle(shift, self.block_shape[1], self.block_shape[0], edgecolor='k', facecolor=next(self.COLORS))]
            else:
                raise AttributeError()

            for i, j in self._positions_of_index(k):
                for patch in patches_of_block:
                    local_shift = np.array(patch.get_xy()) + np.array(((j-k)*self.block_shape[1], i*self.block_shape[0]))
                    patches.append(Rectangle(local_shift, patch.get_width(), patch.get_height(), facecolor=patch.get_facecolor(), alpha=0.5))
            patches.extend(patches_of_block)
        return patches


if __name__ == "__main__":
    # A = BlockSymmetricToeplitzMatrix([np.random.rand(2, 2), np.ones((2, 2)), np.zeros((2, 2))])
    #
    # print(A.nb_blocks)
    # print(A.shape)
    # print(A.full_matrix())
    # A.plot_shape()
    #
    # new_block = lambda: BlockMatrix([
    #     [np.random.rand(2, 2), np.random.rand(2, 2)],
    #     [np.random.rand(2, 2), np.random.rand(2, 2)],
    # ])
    # B = BlockSymmetricToeplitzMatrix([new_block(), new_block()])
    # print(B._t_blocks)
    # print(B.shape)
    # B.plot_shape()

    # new_block = lambda: BlockSymmetricToeplitzMatrix([np.random.rand(2, 2), np.random.rand(2, 2)])
    # C = BlockSymmetricToeplitzMatrix([new_block(), new_block()])
    # print(C._t_blocks)
    # print(C.shape)
    # C.plot_shape()

    random_block = lambda: np.random.rand(1, 1)
    D = BlockSymmetricCirculantMatrix([random_block(), random_block(), random_block(), random_block()])
    print(D._index_grid())
    print(D._positions_of_index(1))
    print(D.shape)
    D.plot_shape()
