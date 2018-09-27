#!/usr/bin/env python
# coding: utf-8

import numpy as np
from capytaine.matrices.block_matrices import BlockMatrix


class BlockSymmetricToeplitzMatrix(BlockMatrix):
    ndim = 2

    def __init__(self, t_blocks):
        # TODO: Check dimensionality
        self.t_blocks = np.asanyarray(t_blocks)

    @property
    def nb_blocks(self):
        return self.t_blocks.shape[0], self.t_blocks.shape[0]

    @property
    def block_shape(self):
        return self.t_blocks[0].shape

    @property
    def shape(self):
        return (sum(block.shape[0] for block in self.t_blocks[:]),
                sum(block.shape[1] for block in self.t_blocks[:]))

    @property
    def T(self):
        transposed_blocks = [block.T for block in self.t_blocks]
        return BlockSymmetricToeplitzMatrix(transposed_blocks)

    @property
    def blocks(self):
        return [[block for block in self.t_blocks[indices]] for indices in self._index_grid()]

    def _index_grid(self):
        n = self.nb_blocks[0]
        return [[abs(i-j) for i in range(n)] for j in range(n)]

    def _positions_of_index(self, k):
        n = self.nb_blocks[0]
        return [(i, j) for i in range(1, n) for j in range(n) if abs(i-j) == k]

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        for k, block in enumerate(self.t_blocks):
            shift = np.array((global_shift[0] + k*self.block_shape[1], global_shift[1]))
            if isinstance(block, BlockMatrix):
                patches_of_block = block._patches(shift)
            elif isinstance(block, np.ndarray):
                patches_of_block = [Rectangle(shift, self.block_shape[1], self.block_shape[0], edgecolor='k', facecolor=next(self.COLORS))]
            else:
                raise AttributeError()

            patches.extend(patches_of_block)
            for i, j in self._positions_of_index(k):
                for patch in patches_of_block:
                    local_shift = np.array(patch.get_xy()) + np.array(((j-k)*self.block_shape[1], i*self.block_shape[0]))
                    patches.append(Rectangle(local_shift, patch.get_width(), patch.get_height(), facecolor=patch.get_facecolor(), alpha=0.5))

        return patches


if __name__ == "__main__":
    A = BlockSymmetricToeplitzMatrix([np.random.rand(2, 2), np.ones((2, 2)), np.zeros((2, 2))])

    print(A.nb_blocks)
    print(A.shape)
    print(A.full_matrix())
    A.plot_shape()

    new_block = lambda: BlockMatrix([
        [np.random.rand(2, 2), np.random.rand(2, 2)],
        [np.random.rand(2, 2), np.random.rand(2, 2)],
    ])
    B = BlockSymmetricToeplitzMatrix([new_block(), new_block()])
    print(B.t_blocks)
    print(B.shape)
    B.plot_shape()

    new_block = lambda: BlockSymmetricToeplitzMatrix([np.random.rand(2, 2), np.random.rand(2, 2)])
    C = BlockSymmetricToeplitzMatrix([new_block(), new_block()])
    print(C.t_blocks)
    print(C.shape)
    C.plot_shape()
