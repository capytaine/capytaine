#!/usr/bin/env python
# coding: utf-8

import numpy as np
from capytaine.matrices.block_matrices import BlockMatrix


class BlockSymmetricToeplitzMatrix:
    ndim = 2

    def __init__(self, t_blocks):
        # TODO: Check dimensionality
        self.t_blocks = np.asanyarray(t_blocks)

    @property
    def nb_blocks(self):
        return (self.t_blocks.shape[0], self.t_blocks.shape[0])

    @property
    def shape(self):
        return (sum(block.shape[0] for block in self.t_blocks[:]),
                sum(block.shape[1] for block in self.t_blocks[:]))

    @property
    def T(self):
        transposed_blocks = [block.T for block in self.t_blocks]
        return BlockSymmetricToeplitzMatrix(transposed_blocks)

    def full_matrix(self):
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in self.t_blocks[indices]]
                       for indices in self._index_grid()]
        return np.block(full_blocks)

    def _index_grid(self):
        n = self.nb_blocks[0]
        return np.fromfunction(lambda i, j: abs(i-j), (n, n), dtype=np.int)


if __name__ == "__main__":
    # print(BlockSymmetricToeplitzMatrix._index_grid(5))

    A = BlockSymmetricToeplitzMatrix([np.random.rand(2, 2), np.ones((2, 2)), np.zeros((2, 2))])

    print(A.nb_blocks)
    print(A.shape)
    print(A.full_matrix())

    B = BlockSymmetricToeplitzMatrix([A, A.T])

    print(B.t_blocks)
    print(B.shape)

