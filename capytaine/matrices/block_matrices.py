#!/usr/bin/env python
# coding: utf-8

import logging

from numbers import Number
from typing import Tuple, List, Callable, Union, Iterable
from itertools import cycle, accumulate, chain, product

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

        self._stored_nb_blocks = self._stored_blocks.shape[:self.ndim]

        # Total shape of the full matrix
        self.shape = self._compute_shape()

        self.dtype = self._stored_blocks[0][0].dtype

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

    def _compute_shape(self):
        return sum(self._stored_block_shapes[0]), sum(self._stored_block_shapes[1])

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

    def _stored_block_positions(self) -> Iterable[List[Tuple[int, int]]]:
        """The position of each blocks in the matrix.

        Example::

            AAAABB
            AAAABB  ->  list(matrix._stored_block_positions) = [[(0,0)], [(0, 4)], [(2, 0)], [(2, 4)]]
            CCCCDD
        """
        x_acc = accumulate([0] + self.block_shapes[0][:-1])
        y_acc = accumulate([0] + self.block_shapes[1][:-1])
        return ([(x, y)] for x, y in product(x_acc, y_acc))

    # TRANSFORMING DATA

    def _apply_unary_op(self, op: Callable) -> 'BlockMatrix':
        """Helper function applying a function recursively on all submatrices."""
        LOG.debug(f"Apply op {op.__name__} to {self}")
        result = [[op(block) for block in line] for line in self._stored_blocks]
        return self.__class__(result, _stored_block_shapes=self._stored_block_shapes, check_dim=False)

    def _apply_binary_op(self, op: Callable, other: 'BlockMatrix') -> 'BlockMatrix':
        """Helper function applying a binary operator recursively on all submatrices."""
        if isinstance(other, self.__class__) and self.nb_blocks == other.nb_blocks:
            LOG.debug(f"Apply op {op.__name__} to {self} and {other}")
            result = [
                [op(block, other_block) for block, other_block in zip(line, other_line)]
                for line, other_line in zip(self._stored_blocks, other._stored_blocks)
            ]
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

    def matvec(self, other):
        """Matrix vector product.
        Named as such to be used as scipy LinearOperator."""
        result = np.zeros(self.shape[0], dtype=other.dtype)
        line_heights = self.block_shapes[0]
        line_positions = list(accumulate(chain([0], line_heights)))
        col_widths = self.block_shapes[1]
        col_positions = list(accumulate(chain([0], col_widths)))
        for line, line_position, line_height in zip(self.all_blocks, line_positions, line_heights):
            line_slice = slice(line_position, line_position+line_height)
            for block, col_position, col_width in zip(line, col_positions, col_widths):
                col_slice = slice(col_position, col_position+col_width)
                result[line_slice] += block @ other[col_slice]
        return result

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
            if other.ndim == 2:
                LOG.debug(f"Multiplication of {self} with a full matrix of shape {other.shape}.")
                # Cut the matrix and recursively call itself to use the code above.
                from capytaine.matrices.builders import cut_matrix
                cut_other = cut_matrix(other, self.block_shapes[1], [other.shape[1]], check_dim=False)
                return (self @ cut_other).full_matrix()

            elif other.ndim == 1:
                LOG.debug(f"Multiplication of {self} with a full vector of size {other.shape}.")
                return self.matvec(other)

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

    def _put_in_full_matrix(self, full_matrix, where=(0, 0)):
        """In place copy the content of the block matrix in a matrix."""
        positions_of_blocks = ([(where[0] + x, where[1] + y) for x, y, in positions_of_block]
                               for positions_of_block in self._stored_block_positions())
        all_blocks_flat = (block for line in self._stored_blocks for block in line)
        for block, positions_of_block in zip(all_blocks_flat, positions_of_blocks):
            if isinstance(block, BlockMatrix):
                position_of_first_appearance = positions_of_block[0]
                block._put_in_full_matrix(full_matrix, where=position_of_first_appearance)
                frame_of_first_appearance = (slice(position_of_first_appearance[0], position_of_first_appearance[0]+block.shape[0]),
                                             slice(position_of_first_appearance[1], position_of_first_appearance[1]+block.shape[1]))
                for position in positions_of_block[1:]:
                    block_frame = (slice(position[0], position[0]+block.shape[0]),
                                   slice(position[1], position[1]+block.shape[1]))
                    full_matrix[block_frame] = full_matrix[frame_of_first_appearance]
            else:
                full_block = block if isinstance(block, np.ndarray) else block.full_matrix()
                for position in positions_of_block:
                    block_frame = (slice(position[0], position[0]+block.shape[0]),
                                   slice(position[1], position[1]+block.shape[1]))
                    full_matrix[block_frame] = full_block
        return full_matrix

    def full_matrix(self) -> np.ndarray:
        """Flatten the block structure and return a full matrix."""
        full_matrix = np.empty(self.shape, dtype=self.dtype)
        self._put_in_full_matrix(full_matrix)
        return full_matrix

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
        blocks_flat_list = (block for line in self.all_blocks for block in line)
        for block_position_in_local_frame, block in zip(self._stored_block_positions(), blocks_flat_list):
            block_position_in_global_frame = (global_frame[0] + block_position_in_local_frame[0][1],
                                              global_frame[1] + block_position_in_local_frame[0][0])
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

