#!/usr/bin/env python
# coding: utf-8

import logging

from numbers import Number
from typing import Tuple, List, Callable, Union
from itertools import cycle

import numpy as np
from matplotlib.patches import Rectangle

LOG = logging.getLogger(__name__)


class BlockMatrix:
    """A (2D) matrix, stored as a set of submatrices (or blocks)."""

    ndim = 2  # Other dimensions have not been implemented.
    display_color = cycle([f'C{i}' for i in range(10)])

    def __init__(self, blocks, _stored_block_shapes=None, check_dim=True):
        assert blocks[0][0].ndim == self.ndim
        self._stored_blocks = np.asarray(blocks)

        if _stored_block_shapes is None:
            self._stored_block_shapes = ([block.shape[0] for block in self._stored_blocks[:, 0]],
                                         [block.shape[1] for block in self._stored_blocks[0, :]])
        else:
            # To avoid going through the whole tree if it is already known.
            self._stored_block_shapes = _stored_block_shapes

        self._stored_shape = (sum(self._stored_block_shapes[0]), sum(self._stored_block_shapes[1]))
        self._stored_nb_blocks = self._stored_blocks.shape[:self.ndim]

        LOG.debug(f"New block matrix: %s", self)

        if check_dim:
            self._check_dimension()

    def __hash__(self):
        # Temporary
        return id(self)

    # ACCESSING DATA

    @property
    def all_blocks(self) -> np.ndarray:
        """The matrix of matrices. For a full block matrix, all the blocks are stored in memory."""
        return self._stored_blocks

    @property
    def shape(self) -> Tuple[int, int]:
        """The total size of the matrix.
        That is the shape of the array if the BlockMatrix was flattened.

        Example::

            AAAABB
            AAAABB  ->  shape = (3, 6)
            AAAABB
        """
        return self._stored_shape

    @property
    def block_shapes(self) -> Tuple[List[int], List[int]]:
        """The shapes of the blocks composing the BlockMatrix.

        Example::

            AAAABB
            AAAABB  ->  block_shapes = ([3], [4, 2])
            AAAABB
        """
        return self._stored_block_shapes

    @property
    def nb_blocks(self) -> Tuple[int, int]:
        """The number of blocks in each directions.

        Example::

            AAAABB
            AAAABB  ->  nb_blocks = (1, 2)
            AAAABB
        """
        return self._stored_nb_blocks

    def _check_dimension(self) -> None:
        """Check that the dimensions of the blocks are consistent."""
        for line in self.all_blocks:
            block_height = line[0].shape[0]
            for block in line[1:]:
                assert block.shape[0] == block_height  # Same height on a given line

        for col in np.moveaxis(self.all_blocks, 1, 0):
            block_width = col[0].shape[1]
            for block in col[1:]:
                assert block.shape[1] == block_width  # Same width on a given column

    @property
    def _block_positions_list(self) -> List[Tuple[int, int]]:
        """The position of each blocks in the matrix.

        Example::

            AAAABB
            AAAABB  ->  _block_positions_list = [(0,0), (0, 4)]
            AAAABB
        """
        positions = []
        current_cursor_position = np.array([0, 0], dtype=np.int)
        for line in self.all_blocks:
            for block in line:
                positions.append(tuple(current_cursor_position))
                current_cursor_position[0] += block.shape[1]  # Shift to the next block on the line
            current_cursor_position[0] = 0  # Return to the beginning of line
            current_cursor_position[1] += line[0].shape[0]  # Shift to next line
        return positions

    # TRANSFORMING DATA

    def _apply_unary_op(self, op: Callable) -> 'BlockMatrix':
        """Helper function applying a function recursively on all submatrices."""
        LOG.debug(f"Apply op {op.__name__} to {self}")
        result = [[op(block) for block in line] for line in self._stored_blocks]
        result = np.asarray(result)
        return self.__class__(result, _stored_block_shapes=self._stored_block_shapes, check_dim=False)

    def _apply_binary_op(self, op: Callable, other: 'BlockMatrix') -> 'BlockMatrix':
        """Helper function applying a binary operator recursively on all submatrices."""
        if isinstance(other, self.__class__) and self.nb_blocks == other.nb_blocks:
            LOG.debug(f"Apply op {op.__name__} to {self} and {other}")
            result = [
                [op(block, other_block) for block, other_block in zip(line, other_line)]
                for line, other_line in zip(self._stored_blocks, other._stored_blocks)
            ]
            result = np.asarray(result)
            return self.__class__(result, _stored_block_shapes=self._stored_block_shapes, check_dim=False)
        else:
            return NotImplemented

    def __add__(self, other: 'BlockMatrix') -> 'BlockMatrix':
        from operator import add
        return self._apply_binary_op(add, other)

    def __radd__(self, other: 'BlockMatrix') -> 'BlockMatrix':
        return self + other

    def __neg__(self) -> 'BlockMatrix':
        from operator import neg
        return self._apply_unary_op(neg)

    def __sub__(self, other: 'BlockMatrix') -> 'BlockMatrix':
        return self + (-other)

    def __rsub__(self, other: 'BlockMatrix') -> 'BlockMatrix':
        return other + (-self)

    def __mul__(self, other: Union['BlockMatrix', Number]) -> 'BlockMatrix':
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: other*x)
        else:
            from operator import mul
            return self._apply_binary_op(mul, other)

    def __rmul__(self, other: Union['BlockMatrix', Number]) -> 'BlockMatrix':
        return self * other

    def __truediv__(self, other: Union['BlockMatrix', Number]) -> 'BlockMatrix':
        from numbers import Number
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: x/other)
        else:
            from operator import truediv
            return self._apply_binary_op(truediv, other)

    def __rtruediv__(self, other: Union['BlockMatrix', Number]) -> 'BlockMatrix':
        from numbers import Number
        if isinstance(other, Number):
            return self._apply_unary_op(lambda x: other/x)
        else:
            return self._apply_binary_op(lambda x, y: y/x, other)

    def __matmul__(self, other: Union['BlockMatrix', np.ndarray]) -> Union['BlockMatrix', np.ndarray]:
        if isinstance(other, BlockMatrix) and self.block_shapes[1] == other.block_shapes[0]:
            LOG.debug(f"Multiplication of %s with %s", self, other)
            own_blocks = self.all_blocks
            other_blocks = np.moveaxis(other.all_blocks, 1, 0)
            new_matrix = []
            for own_line in own_blocks:
                new_line = []
                for other_col in other_blocks:
                    new_line.append(sum(own_block @ other_block for own_block, other_block in zip(own_line, other_col)))
                new_matrix.append(new_line)
            return BlockMatrix(new_matrix, check_dim=False)

        elif isinstance(other, np.ndarray) and self.shape[1] == other.shape[0]:
            LOG.debug(f"Multiplication of %s with a vector.", self)
            if other.ndim == 2:
                from capytaine.matrices.builders import cut_matrix
                cut_other = cut_matrix(other, self.block_shapes[1], [other.shape[1]], check_dim=False)
                return (self @ cut_other).full_matrix()
            elif other.ndim == 1:
                other = other.reshape((other.shape[0], 1))
                return (self @ other).flatten()
            else:
                return NotImplemented

        else:
            return NotImplemented

    def astype(self, dtype: np.dtype) -> 'BlockMatrix[dtype]':
        return self._apply_unary_op(lambda x: x.astype(dtype))

    @property
    def T(self) -> 'BlockMatrix':
        """Transposed matrix."""
        transposed_blocks = np.array([[block.T for block in line] for line in self.all_blocks])
        return BlockMatrix(transposed_blocks.T, check_dim=False)

    def full_matrix(self) -> np.ndarray:
        """Flatten the block structure and return a full matrix."""
        full_blocks = [[block.full_matrix() if not isinstance(block, np.ndarray) else block
                        for block in line]
                       for line in self.all_blocks]
        return np.block(full_blocks)

    # COMPARISON AND REDUCTION

    def __eq__(self, other: 'BlockMatrix') -> 'BlockMatrix[bool]':
        from operator import eq
        return self._apply_binary_op(eq, other)

    def __invert__(self) -> 'BlockMatrix':
        """Boolean not (~)"""
        from operator import invert
        return self._apply_unary_op(invert)

    def __ne__(self, other: 'BlockMatrix') -> 'BlockMatrix[bool]':
        return ~(self == other)

    def all(self) -> bool:
        for line in self._stored_blocks:
            for block in line:
                if not block.all():
                    return False
        return True

    def any(self) -> bool:
        for line in self._stored_blocks:
            for block in line:
                if block.any():
                    return True
        return False

    def min(self) -> Number:
        return min(block.min() for line in self._stored_blocks for block in line)

    def max(self) -> Number:
        return max(block.max() for line in self._stored_blocks for block in line)

    # DISPLAYING DATA

    def __repr__(self):
        return f"{self.__class__.__name__}(nb_blocks={self.nb_blocks}, shape={self.shape})"

    def _patches(self, global_frame: Tuple[int, int]) -> List[Rectangle]:
        """Helper function for displaying the shape of the matrix.
        Recursively returns a list of rectangles representing the sub-blocks of the matrix.

        Parameters
        ----------
        global_frame: tuple of ints
            coordinates of the origin in the top left corner.
        """
        patches = []
        for block_position_in_local_frame, block in zip(self._block_positions_list, [block for line in self.all_blocks for block in line]):
            block_position_in_global_frame = (global_frame[0] + block_position_in_local_frame[0],
                                              global_frame[1] + block_position_in_local_frame[1])
            if isinstance(block, BlockMatrix):
                patches.extend(block._patches(block_position_in_global_frame))
            elif isinstance(block, np.ndarray):
                patches.append(Rectangle(block_position_in_global_frame, block.shape[1], block.shape[0],
                                         edgecolor='k', facecolor=next(self.display_color)))
            else:
                raise AttributeError()
        return patches

    def plot_shape(self):
        """Plot the structure of the matrix using matplotlib."""
        import matplotlib.pyplot as plt
        plt.figure()
        for patch in self._patches((0, 0)):
            plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.xlim(0, self.shape[1])
        plt.ylim(0, self.shape[0])
        plt.gca().invert_yaxis()
        plt.show()

