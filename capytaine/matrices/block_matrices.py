#!/usr/bin/env python
# coding: utf-8

from itertools import cycle

import numpy as np


class BlockMatrix:
    """A (2D) matrix, stored as a set of submatrices (or blocks)."""

    ndim = 2  # Other dimensions have not been implemented.
    display_color = cycle([f'C{i}' for i in range(10)])

    def __init__(self, blocks):
        assert blocks[0][0].ndim == self.ndim
        self._stored_blocks = np.asanyarray(blocks)
        self._check_dimension()

    # ACCESSING DATA

    @property
    def all_blocks(self):
        return self._stored_blocks

    @property
    def all_blocks_list(self):
        return [block for line in self.all_blocks for block in line]

    @property
    def _stored_blocks_flat(self):
        return np.array([block for line in self._stored_blocks for block in line])

    def _check_dimension(self) -> None:
        for line in self.all_blocks:
            for block in line:
                assert block.shape[0] == line[0].shape[0]  # Same height on a given line

        for col in np.moveaxis(self.all_blocks, 1, 0):
            for block in col:
                assert block.shape[1] == col[0].shape[1]  # Same width on a given column

    @property
    def shape(self):
        return (sum(block.shape[0] for block in self.all_blocks[:, 0]),
                sum(block.shape[1] for block in self.all_blocks[0, :]))

    @property
    def nb_blocks(self):
        return self.all_blocks.shape[:self.ndim]

    @property
    def _block_positions_list(self):
        positions = []
        position = np.array([0, 0], dtype=np.int)
        for line in self.all_blocks:
            for block in line:
                positions.append(position.copy())
                position[0] += block.shape[1]  # Shift to the next block on the line
            position[0] = 0  # Return to the beginning of line
            position[1] += line[0].shape[0]  # Shift to next line
        return positions

    # DISPLAYING DATA

    def __repr__(self):
        return f"{self.__class__.__name__}(nb_blocks={self.nb_blocks}, shape={self.shape})"

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        patches = []
        for position, block in zip(self._block_positions_list, self.all_blocks_list):
            if isinstance(block, BlockMatrix):
                patches.extend(block._patches(global_shift + position))
            elif isinstance(block, np.ndarray):
                patches.append(Rectangle(global_shift + position, block.shape[1], block.shape[0],
                                         edgecolor='k', facecolor=next(self.display_color)))
            else:
                raise AttributeError()
        return patches

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

    # TRANSFORMING DATA

    @property
    def T(self):
        transposed_blocks = [[block.T for block in line] for line in self.all_blocks]
        return BlockMatrix(transposed_blocks)

    def full_matrix(self):
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in line]
                       for line in self.all_blocks]
        return np.block(full_blocks)


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
