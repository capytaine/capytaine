#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


class BlockMatrix:
    ndim = 2  # Other dimensions have not been implemented.

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

    def __repr__(self):
        return f"{self.__class__.__name__}(nb_blocks={self.nb_blocks}, shape={self.shape})"

    @property
    def T(self):
        transposed_blocks = [[block.T for block in line] for line in self.blocks]
        return BlockMatrix(transposed_blocks)

    def full_matrix(self):
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in line]
                       for line in self.blocks]
        return np.block(full_blocks)

    def plot_shape(self):
        fig = plt.figure()
        coll = PatchCollection(self._patches(np.zeros((2,))))
        plt.gca().add_collection(coll)
        # plt.axis('off')
        plt.show()

    def _patches(self, global_shift):
        patches = []
        shift = global_shift.copy()
        for line in self.blocks:
            shift[0] = global_shift[0]
            for block in line:
                if isinstance(block, np.ndarray):
                    patches.append(Rectangle(shift, block.shape[0], block.shape[1]))
                    shift[0] += block.shape[0]
            shift[0] += block.shape[1]
        return patches


def zeros_like(A):
    if isinstance(A, BlockMatrix):
        I = []
        for i in range(A.nb_blocks[0]):
            line = []
            for j in range(A.nb_blocks[1]):
                line.append(zeros_like(A.blocks[i][j]))
            I.append(line)
        return BlockMatrix(I)
    elif isinstance(A, np.ndarray):
        return np.zeros_like(A)


def identity_like(A):
    if isinstance(A, BlockMatrix):
        I = []
        for i in range(A.nb_blocks[0]):
            line = []
            for j in range(A.nb_blocks[1]):
                if i == j:
                    line.append(identity_like(A.blocks[i][j]))
                else:
                    line.append(zeros_like(A.blocks[i][j]))
            I.append(line)
        return BlockMatrix(I)
    elif isinstance(A, np.ndarray):
        return np.eye(A.shape[0], A.shape[1])


if __name__ == "__main__":

    A = BlockMatrix([
        [np.random.rand(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.random.rand(2, 2)]
    ])
    print(repr(A))
    print(A.full_matrix())
    A.plot_shape()


    # I = identity_like(A)
    # print(repr(I))
    # print(I.full_matrix())

    # B = BlockMatrix([[A, np.zeros((4, 1))]])
    # print(repr(B))
    # print(B.full_matrix())

    # O = zeros_like(B)
    # print(repr(O))
    # print(O.full_matrix())

    # C = BlockMatrix([
    #     [np.random.rand(3, 3), np.random.rand(3, 1)],
    #     [np.random.rand(1, 3), np.random.rand(1, 1)]
    # ])
    # print(repr(C))
    # print(C.full_matrix())

