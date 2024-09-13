"""This module implements block matrices to be used in Hierarchical Toeplitz matrices.

It takes inspiration from the following works:

* `openHmx module from Gypsilab by Matthieu Aussal (GPL licensed) <https://github.com/matthieuaussal/gypsilab>`_
* `HierarchicalMatrices by Markus Neumann (GPL licensed) <https://github.com/maekke97/HierarchicalMatrices>`_
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

from numbers import Number
from typing import Tuple, List, Callable, Union, Iterable
from itertools import cycle, accumulate, chain, product
from collections.abc import Iterator

import numpy as np

from capytaine.matrices.low_rank import LowRankMatrix
from capytaine.tools.optional_imports import import_optional_dependency

LOG = logging.getLogger(__name__)


class BlockMatrix:
    """A (2D) matrix, stored as a set of submatrices (or blocks).

    Parameters
    ----------
    blocks: list of list of matrices
        The blocks of the block matrix.
    check: bool, optional
        Should the code perform sanity checks on the inputs? (default: True)

    Attributes
    ----------
    shape: pair of ints
        shape of the full matrix
    nb_blocks: pair of ints
        number of blocks in each directions.
        Example::

            AAAABB
            AAAABB  ->  nb_blocks = (1, 2)
            AAAABB
    """

    ndim = 2  # Other dimensions have not been implemented.

    def __init__(self, blocks, _stored_block_shapes=None, check=True):

        self._stored_nb_blocks = (len(blocks), len(blocks[0]))

        self._stored_blocks = np.empty(self._stored_nb_blocks, dtype=object)
        for i in range(len(blocks)):
            for j in range(len(blocks[i])):
                self._stored_blocks[i, j] = blocks[i][j]

        if _stored_block_shapes is None:
            self._stored_block_shapes = ([block.shape[0] for block in self._stored_blocks[:, 0]],
                                         [block.shape[1] for block in self._stored_blocks[0, :]])
        else:
            # To avoid going down the tree if it is already known.
            self._stored_block_shapes = _stored_block_shapes

        # Block shape of the full matrix
        self.nb_blocks = self._compute_nb_blocks()

        # Total shape of the full matrix
        self.shape = self._compute_shape()

        LOG.debug(f"New block matrix: %s", self)

        if check:
            assert self._check_dimensions_of_blocks()
            assert self._check_dtype()

    def _compute_shape(self):
        # In a dedicated routine because it will be overloaded by subclasses.
        return sum(self._stored_block_shapes[0]), sum(self._stored_block_shapes[1])

    def _compute_nb_blocks(self):
        return self._stored_nb_blocks

    def _check_dimensions_of_blocks(self) -> bool:
        """Check that the dimensions of the blocks are consistent."""
        if not all(block.ndim == self.ndim for line in self._stored_blocks for block in line):
            return False

        for line in self.all_blocks:
            block_height = line[0].shape[0]
            for block in line[1:]:
                if not block.shape[0] == block_height:  # Same height on a given line
                    return False

        for col in np.moveaxis(self.all_blocks, 1, 0):
            block_width = col[0].shape[1]
            for block in col[1:]:
                if not block.shape[1] == block_width:  # Same width on a given column
                    return False
        return True

    def _check_dtype(self) -> bool:
        """Check that the type of the blocks are consistent."""
        for line in self._stored_blocks:
            for block in line:
                if block.dtype != self.dtype:
                    return False
        return True

    # ACCESSING DATA

    @property
    def dtype(self):
        """The type of data of all of the subblocks."""
        try:
            return self._stored_blocks[0][0].dtype
        except AttributeError:
            return None

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

    def _stored_block_positions(self, global_frame=(0, 0)) -> Iterable[List[Tuple[int, int]]]:
        """The position of each blocks in the matrix as a generator.
        The list is used by subclasses where the same block may appear several times in different positions.

        Parameters
        ----------
        global_frame: Tuple[int], optional
            the coordinate of the top right corner. Default: 0, 0.

        Example::

            AAAABB
            AAAABB  ->  list(matrix._stored_block_positions) = [[(0,0)], [(0, 4)], [(2, 0)], [(2, 4)]]
            CCCCDD
        """
        x_acc = accumulate([0] + self.block_shapes[0][:-1])
        y_acc = accumulate([0] + self.block_shapes[1][:-1])
        return ([(global_frame[0] + x, global_frame[1] + y)] for x, y in product(x_acc, y_acc))

    def _put_in_full_matrix(self, full_matrix: np.ndarray, where=(0, 0)) -> None:
        """Copy the content of the block matrix in a matrix, which is modified in place."""
        all_blocks_in_flat_iterator = (block for line in self._stored_blocks for block in line)
        positions_of_all_blocks = self._stored_block_positions(global_frame=where)
        for block, positions_of_the_block in zip(all_blocks_in_flat_iterator, positions_of_all_blocks):
            if isinstance(block, BlockMatrix):
                position_of_first_appearance = positions_of_the_block[0]
                frame_of_first_appearance = (slice(position_of_first_appearance[0], position_of_first_appearance[0]+block.shape[0]),
                                             slice(position_of_first_appearance[1], position_of_first_appearance[1]+block.shape[1]))
                block._put_in_full_matrix(full_matrix, where=position_of_first_appearance)

                for position in positions_of_the_block[1:]:  # For the other appearances, only copy the first appearance
                    block_frame = (slice(position[0], position[0]+block.shape[0]),
                                   slice(position[1], position[1]+block.shape[1]))
                    full_matrix[block_frame] = full_matrix[frame_of_first_appearance]

            else:
                full_block = block if isinstance(block, np.ndarray) else block.full_matrix()
                for position in positions_of_the_block:
                    block_frame = (slice(position[0], position[0]+block.shape[0]),
                                   slice(position[1], position[1]+block.shape[1]))
                    full_matrix[block_frame] = full_block

    def full_matrix(self, dtype=None) -> np.ndarray:
        """Flatten the block structure and return a full matrix."""
        if dtype is None: dtype = self.dtype
        full_matrix = np.empty(self.shape, dtype=dtype)
        self._put_in_full_matrix(full_matrix)
        return full_matrix

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise ValueError("Making an ndarray out of a BlockMatrix requires copy")
        return self.full_matrix(dtype=dtype)

    def no_toeplitz(self):
        """Recursively replace the block toeplitz matrices by usual block matrices.
        WARNING: the block matrices may still contain several references to the same block."""
        blocks = [[block.no_toeplitz() if isinstance(block, BlockMatrix) else block for block in line] for line in self.all_blocks]
        return BlockMatrix(blocks)

    def __deepcopy__(self, memo):
        from copy import deepcopy
        blocks = [[deepcopy(block) for block in line] for line in self._stored_blocks]
        # The call to deepcopy does not use the memo on purpose:
        # the goal is to replace references to the same block by references to different copies of the block.
        return self.__class__(blocks)

    @property
    def stored_data_size(self):
        """Return the number of entries actually stored in memory."""
        size = 0
        for line in self._stored_blocks:
            for block in line:
                if isinstance(block, np.ndarray):
                    size += np.prod(block.shape)
                else:
                    size += block.stored_data_size
        return size

    @property
    def density(self):
        return self.stored_data_size/np.prod(self.shape)

    @property
    def sparcity(self):
        return 1 - self.density

    def __hash__(self):
        # Temporary
        return id(self)

    # TRANSFORMING DATA

    def _apply_unary_op(self, op: Callable) -> 'BlockMatrix':
        """Helper function applying a function recursively on all submatrices."""
        LOG.debug(f"Apply op {op.__name__} to {self}")
        result = [[op(block) for block in line] for line in self._stored_blocks]
        return self.__class__(result, _stored_block_shapes=self._stored_block_shapes, check=False)

    def _apply_binary_op(self, op: Callable, other: 'BlockMatrix') -> 'BlockMatrix':
        """Helper function applying a binary operator recursively on all submatrices."""
        if isinstance(other, self.__class__) and self.nb_blocks == other.nb_blocks:
            LOG.debug(f"Apply op {op.__name__} to {self} and {other}")
            result = [
                [op(block, other_block) for block, other_block in zip(line, other_line)]
                for line, other_line in zip(self._stored_blocks, other._stored_blocks)
            ]
            return self.__class__(result, _stored_block_shapes=self._stored_block_shapes, check=False)
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
        from operator import sub
        return self._apply_binary_op(sub, other)

    def __rsub__(self, other: 'BlockMatrix') -> 'BlockMatrix':
        from operator import sub
        return other._apply_binary_op(sub, self)

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
        LOG.debug(f"Multiplication of {self} with a full vector of size {other.shape}.")
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

    def rmatvec(self, other):
        """Vector matrix product.
        Named as such to be used as scipy LinearOperator."""
        LOG.debug(f"Multiplication of a full vector of size {other.shape} with {self}.")
        result = np.zeros(self.shape[1], dtype=other.dtype)
        line_heights = self.block_shapes[0]
        line_positions = list(accumulate(chain([0], line_heights)))
        col_widths = self.block_shapes[1]
        col_positions = list(accumulate(chain([0], col_widths)))
        for col, col_position, col_width in zip(self.all_blocks.T, col_positions, col_widths):
            col_slice = slice(col_position, col_position+col_width)
            for block, line_position, line_height in zip(col, line_positions, line_heights):
                line_slice = slice(line_position, line_position+line_height)
                if isinstance(block, BlockMatrix):
                    result[col_slice] += block.rmatvec(other[line_slice])
                else:
                    result[col_slice] += other[line_slice] @ block
        return result

    def matmat(self, other):
        """Matrix-matrix product."""
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
            return BlockMatrix(new_matrix, check=False)

        elif isinstance(other, np.ndarray) and self.shape[1] == other.shape[0]:
            LOG.debug(f"Multiplication of {self} with a full matrix of shape {other.shape}.")
            # Cut the matrix and recursively call itself to use the code above.
            from capytaine.matrices.builders import cut_matrix
            cut_other = cut_matrix(other, self.block_shapes[1], [other.shape[1]], check=False)
            return (self @ cut_other).full_matrix()

    def __matmul__(self, other: Union['BlockMatrix', np.ndarray]) -> Union['BlockMatrix', np.ndarray]:
        if not (isinstance(other, BlockMatrix) or isinstance(other, np.ndarray)):
            return NotImplemented
        elif other.ndim == 2:  # Other is a matrix
            if other.shape[1] == 1:  # Actually a column vector
                return self.matvec(other.flatten())
            else:
                return self.matmat(other)
        elif other.ndim == 1:  # Other is a vector
            return self.matvec(other)
        else:
            return NotImplemented

    def __rmatmul__(self, other: Union['BlockMatrix', np.ndarray]) -> Union['BlockMatrix', np.ndarray]:
        if not (isinstance(other, BlockMatrix) or isinstance(other, np.ndarray)):
            return NotImplemented
        elif other.ndim == 2:  # Other is a matrix
            if other.shape[1] == 1:  # Actually a column vector
                return self.rmatvec(other.flatten())
            else:
                return NotImplemented
        elif other.ndim == 1:  # Other is a vector
            return self.rmatvec(other)
        else:
            return NotImplemented

    def astype(self, dtype: np.dtype) -> 'BlockMatrix':
        return self._apply_unary_op(lambda x: x.astype(dtype))

    def fft_of_list(*block_matrices, check=True):
        """Compute the fft of a list of block matrices of the same type and shape.
        The output is a list of block matrices of the same shape as the input ones.
        The fft is computed element-wise, so the block structure does not cause any mathematical difficulty.
        Returns an array of BlockMatrices.
        """
        class_of_matrices = type(block_matrices[0])
        nb_blocks = block_matrices[0]._stored_nb_blocks

        LOG.debug(f"FFT of {len(block_matrices)} {class_of_matrices.__name__} (stored blocks = {nb_blocks})")

        if check:
            # Check the validity of the shapes of the matrices given as input
            shape = block_matrices[0].shape
            assert [nb_blocks == matrix._stored_nb_blocks for matrix in block_matrices[1:]]
            assert [shape == matrix.shape for matrix in block_matrices[1:]]
            assert [class_of_matrices == type(matrix) for matrix in block_matrices[1:]]

        # Initialize a vector of block matrices without values in the blocks.
        result = np.empty(len(block_matrices), dtype=object)
        for i in range(len(block_matrices)):
            result[i] = class_of_matrices(np.empty(nb_blocks, dtype=object),
                                          _stored_block_shapes=block_matrices[0]._stored_block_shapes,
                                          check=False)

        for i_block, j_block in product(range(nb_blocks[0]), range(nb_blocks[1])):
            list_of_i_j_blocks = [block_matrices[i_matrix]._stored_blocks[i_block, j_block]
                                  for i_matrix in range(len(block_matrices))]

            if any(isinstance(block, np.ndarray) or isinstance(block, LowRankMatrix) for block in list_of_i_j_blocks):
                list_of_i_j_blocks = [block if isinstance(block, np.ndarray) else block.full_matrix() for block in list_of_i_j_blocks]
                fft_of_blocks = np.fft.fft(list_of_i_j_blocks, axis=0)
            else:
                fft_of_blocks = BlockMatrix.fft_of_list(*list_of_i_j_blocks, check=False)

            for matrix, computed_block in zip(result, fft_of_blocks):
                matrix._stored_blocks[i_block, j_block] = computed_block

        return result

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

    @property
    def str_shape(self):
        blocks_str = []
        for line in self.all_blocks:
            for block in line:
                if isinstance(block, BlockMatrix):
                    blocks_str.append(block.str_shape)
                elif isinstance(block, np.ndarray) or isinstance(block, LowRankMatrix):
                    blocks_str.append("{}×{}".format(*block.shape))
                else:
                    blocks_str.append("?×?")

        if len(set(blocks_str)) == 1:
            return "{}×{}×[".format(*self.nb_blocks) + blocks_str[0] + "]"
        else:
            blocks_str = np.array(blocks_str).reshape(self.nb_blocks).tolist()
            return str(blocks_str).replace("'", "")

    def __str__(self):
        if not hasattr(self, '_str'):
            args = [self.str_shape]
            if self.dtype not in [np.float64, float]:
                args.append(f"dtype={self.dtype}")
            self._str = f"{self.__class__.__name__}(" + ", ".join(args) + ")"
        return self._str

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    display_color = cycle([f'C{i}' for i in range(10)])

    def _patches(self,
                 global_frame: Union[Tuple[int, int], np.ndarray]
                 ):
        """Helper function for displaying the shape of the matrix.
        Recursively returns a list of rectangles representing the sub-blocks of the matrix.

        Uses BlockMatrix.display_color to assign color to the blocks.
        By default, it cycles through matplotlib default colors.
        But if display_color is redefined as a callable, it is called with the block as argument.

        Parameters
        ----------
        global_frame: tuple of ints
            coordinates of the origin in the top left corner.

        Returns
        -------
        list of matplotlib.patches.Rectangle
        """
        matplotlib_patches = import_optional_dependency("matplotlib.patches", "matplotlib")
        Rectangle = matplotlib_patches.Rectangle

        all_blocks_in_flat_iterator = (block for line in self._stored_blocks for block in line)
        positions_of_all_blocks = self._stored_block_positions(global_frame=global_frame)
        patches = []
        for block, positions_of_the_block in zip(all_blocks_in_flat_iterator, positions_of_all_blocks):
            position_of_first_appearance = positions_of_the_block[0]
            # Exchange coordinates: row index i -> y, column index j -> x
            position_of_first_appearance = np.array((position_of_first_appearance[1], position_of_first_appearance[0]))

            if isinstance(block, BlockMatrix):
                patches_of_this_block = block._patches(np.array((position_of_first_appearance[1], position_of_first_appearance[0])))
            elif isinstance(block, np.ndarray):

                if isinstance(self.display_color, Iterator):
                    color = next(self.display_color)
                elif callable(self.display_color):
                    color = self.display_color(block)
                else:
                    color = np.random.rand(3)

                patches_of_this_block = [Rectangle(position_of_first_appearance,
                                                   block.shape[1], block.shape[0],
                                                   edgecolor='k', facecolor=color)]
            elif isinstance(block, LowRankMatrix):

                if isinstance(self.display_color, Iterator):
                    color = next(self.display_color)
                elif callable(self.display_color):
                    color = self.display_color(block)
                else:
                    color = np.random.rand(3)

                patches_of_this_block = [
                    # Left block
                    Rectangle(position_of_first_appearance,
                              block.left_matrix.shape[1], block.left_matrix.shape[0],
                              edgecolor='k', facecolor=color),
                    # Top block
                    Rectangle(position_of_first_appearance,
                              block.right_matrix.shape[1], block.right_matrix.shape[0],
                              edgecolor='k', facecolor=color),
                    # Rest of the matrix
                    Rectangle(position_of_first_appearance,
                              block.right_matrix.shape[1], block.left_matrix.shape[0],
                              facecolor=color, alpha=0.2),
                ]
            else:
                raise NotImplementedError()

            patches.extend(patches_of_this_block)

            # For the other appearances, copy the patches of the first appearance
            for block_position in positions_of_the_block[1:]:
                block_position = np.array((block_position[1], block_position[0]))
                for patch in patches_of_this_block:  # A block can be made of several patches.
                    shift = block_position - position_of_first_appearance
                    patch_position = np.array(patch.get_xy()) + shift
                    patches.append(Rectangle(patch_position, patch.get_width(), patch.get_height(),
                                             facecolor=patch.get_facecolor(), alpha=0.2))

        return patches

    def plot_shape(self):
        """Plot the structure of the matrix using matplotlib."""
        matplotlib = import_optional_dependency("matplotlib")
        plt = matplotlib.pyplot

        plt.figure()
        for patch in self._patches((0, 0)):
            plt.gca().add_patch(patch)
        plt.axis('equal')
        plt.xlim(0, self.shape[1])
        plt.ylim(0, self.shape[0])
        plt.gca().invert_yaxis()
        # plt.show()


    def access_block_by_path(self, path):
        """
        Access a diagonal block in a block matrix from the path of the
        corresponding leaf
        """
        this_block = self
        for index in path:
            this_block = this_block.all_blocks[index, index]
        return this_block
