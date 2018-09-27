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
    def block_shapes(self):
        return ([block.shape[0] for block in self.all_blocks[:, 0]],
                [block.shape[1] for block in self.all_blocks[0, :]])

    @property
    def shape(self):
        x_shape, y_shape = self.block_shapes
        return sum(x_shape), sum(y_shape)

    @property
    def nb_blocks(self):
        return self.all_blocks.shape[:self.ndim]

    @property
    def _nb_stored_blocks(self):
        return self._stored_blocks.shape[:self.ndim]

    @property
    def _block_positions_list(self):
        positions = []
        position = np.array([0, 0], dtype=np.int)
        for line in self.all_blocks:
            for block in line:
                positions.append(tuple(position))
                position[0] += block.shape[1]  # Shift to the next block on the line
            position[0] = 0  # Return to the beginning of line
            position[1] += line[0].shape[0]  # Shift to next line
        return positions

    # TRANSFORMING DATA

    def _apply_unary_op(self, op):
        result = np.array([op(block) for block in self._stored_blocks_flat])
        return self.__class__(result.reshape(self._nb_stored_blocks + result.shape[1:]))

    def _apply_binary_op(self, op, other):
        if isinstance(other, self.__class__) and self.nb_blocks == other.nb_blocks:
            result = [op(block, other_block) for block, other_block in zip(self._stored_blocks_flat, other._stored_blocks_flat)]
            result = np.array(result)
            return self.__class__(result.reshape(self._nb_stored_blocks + result.shape[1:]))
        else:
            return NotImplemented

    def __add__(self, other):
        from operator import add
        return self._apply_binary_op(add, other)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        from operator import neg
        return self._apply_unary_op(neg)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        from numbers import Number
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: other*x)
        else:
            from operator import mul
            return self._apply_binary_op(mul, other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        from numbers import Number
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: x/other)
        else:
            from operator import truediv
            return self._apply_binary_op(truediv, other)

    def __rtruediv__(self, other):
        from numbers import Number
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: other/x)
        else:
            return self._apply_binary_op(lambda x, y: y/x, other)

    def __matmul__(self, other):
        if isinstance(other, BlockMatrix) and self.block_shapes[1] == other.block_shapes[0]:
            own_blocks = self.all_blocks
            other_blocks = np.moveaxis(other.all_blocks, 1, 0)
            new_matrix = []
            for own_line in own_blocks:
                new_line = []
                for other_col in other_blocks:
                    new_line.append(sum(own_block @ other_block for own_block, other_block in zip(own_line, other_col)))
                new_matrix.append(new_line)
            return BlockMatrix(new_matrix)

        elif isinstance(other, np.ndarray) and self.shape[1] == other.shape[0]:
            if other.ndim == 2:
                from capytaine.matrices.builders import cut_matrix
                cut_other = cut_matrix(other, self.block_shapes[1], [other.shape[1]])
                return (self @ cut_other).full_matrix()
            elif other.ndim == 1:
                other = other.reshape((other.shape[0], 1))
                return (self @ other).flatten()
            else:
                return NotImplemented

        else:
            return NotImplemented


    def astype(self, dtype):
        return self._apply_unary_op(lambda x: x.astype(dtype))

    @property
    def T(self):
        transposed_blocks = np.array([[block.T for block in line] for line in self.all_blocks])
        return BlockMatrix(transposed_blocks.T)

    def full_matrix(self):
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in line]
                       for line in self.all_blocks]
        return np.block(full_blocks)

    # COMPARISON AND REDUCTION

    def __eq__(self, other):
        from operator import eq
        return self._apply_binary_op(eq, other)

    def __invert__(self):
        from operator import invert
        return self._apply_unary_op(invert)

    def __ne__(self, other):
        return ~(self == other)

    def all(self):
        for block in self._stored_blocks_flat:
            if not block.all():
                return False
            return True

    def any(self):
        for block in self._stored_blocks_flat:
            if block.any():
                return True
            return False

    def min(self):
        return min([block.min() for block in self._stored_blocks_flat])

    def max(self):
        return max([block.max() for block in self._stored_blocks_flat])

    # DISPLAYING DATA

    def __repr__(self):
        return f"{self.__class__.__name__}(nb_blocks={self.nb_blocks}, shape={self.shape})"

    def _patches(self, global_shift):
        from matplotlib.patches import Rectangle
        global_shift = np.asarray(global_shift)
        patches = []
        for position, block in zip(self._block_positions_list, [block for line in self.all_blocks for block in line]):
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

