#!/usr/bin/env python
# coding: utf-8

import logging
import numpy as np


class BlockMatrix:
    ndim = 2

    def __init__(self, blocks):
        # TODO: Check dimensionality
        self.blocks = np.asanyarray(blocks)

    @property
    def nb_blocks(self):
        return self.blocks.shape[:self.ndim]

    @property
    def shape(self):
        return (sum(block.shape[0] for block in self.blocks[:, 0]),
                sum(block.shape[1] for block in self.blocks[0, :]))

    @property
    def T(self):
        transposed_blocks = [[block.T for block in line] for line in self.blocks]
        return BlockMatrix(transposed_blocks)

    def full_matrix(self):
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in line]
                       for line in self.blocks]
        return np.block(full_blocks)


if __name__ == "__main__":

    A = BlockMatrix([
        [np.random.rand(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.random.rand(2, 2)]
    ])
    print(A.blocks.shape)
    print(A.nb_blocks)
    print(A.shape)
    print(A.full_matrix())

    B = BlockMatrix([[A, np.zeros((4, 1))]])
    print(B.blocks.shape)
    print(B.nb_blocks)
    print(B.shape)
    print(B.full_matrix())

    C = BlockMatrix([
        [np.random.rand(3, 3), np.random.rand(3, 1)],
        [np.random.rand(1, 3), np.random.rand(1, 1)]
    ])
    print(C.blocks)
    print(C.nb_blocks)
    print(C.shape)
    print(C.full_matrix())

