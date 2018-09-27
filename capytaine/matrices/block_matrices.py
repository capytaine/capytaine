#!/usr/bin/env python
# coding: utf-8

from itertools import cycle

import numpy as np


class BlockMatrix:
    ndim = 2  # Other dimensions have not been implemented.
    COLORS = cycle([f'C{i}' for i in range(10)])

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
        import matplotlib.pyplot as plt
        plt.figure()
        for patch in self._patches(np.zeros((2,))):
            plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.xlim(0, self.shape[1])
        plt.ylim(0, self.shape[0])
        plt.gca().invert_yaxis()
        plt.show()

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        shift = global_shift.copy()
        for line in self.blocks:
            shift[0] = global_shift[0]
            for block in line:
                if isinstance(block, BlockMatrix):
                    patches.extend(block._patches(shift))
                elif isinstance(block, np.ndarray):
                    patches.append(Rectangle(shift, block.shape[1], block.shape[0],
                                             edgecolor='k', facecolor=next(self.COLORS)))
                else:
                    raise AttributeError()
                shift[0] += block.shape[1]
            shift[1] += line[0].shape[0]
        return patches


if __name__ == "__main__":

    A = BlockMatrix([
        [np.random.rand(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.random.rand(2, 2)]
    ])
    print(repr(A))
    print(A.full_matrix())
    A.plot_shape()

    B = BlockMatrix([[A, np.zeros((4, 1))]])
    print(repr(B))
    print(B.full_matrix())
    B.plot_shape()

    C = BlockMatrix([
        [np.random.rand(3, 3), np.random.rand(3, 1)],
        [np.random.rand(1, 3), np.random.rand(1, 1)]
    ])
    print(repr(C))
    print(C.full_matrix())
    C.plot_shape()
