#!/usr/bin/env python
# coding: utf-8
"""This module contains some helpful functions to create block matrices."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from itertools import accumulate

import numpy as np

from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.low_rank import LowRankMatrix

LOG = logging.getLogger(__name__)


def cut_matrix(full_matrix, x_shapes, y_shapes, check=False):
    """Transform a numpy array into a block matrix of numpy arrays.

    Parameters
    ----------
    full_matrix: numpy array
        The matrix to split into blocks.
    x_shapes: sequence of int
        The columns at which to split the blocks.
    y_shapes: sequence of int
        The lines at which to split the blocks.
    check: bool, optional
        Check to dimensions and type of the matrix after creation (default: False).

    Return
    ------
    BlockMatrix
        The same matrix as the input one but in block form.
    """
    new_block_matrix = []
    for i, di in zip(accumulate([0] + x_shapes[:-1]), x_shapes):
        line = []
        for j, dj in zip(accumulate([0] + x_shapes[:-1]), y_shapes):
            line.append(full_matrix[i:i+di, j:j+dj])
        new_block_matrix.append(line)
    return BlockMatrix(new_block_matrix, check=check)


def random_block_matrix(x_shapes, y_shapes):
    """A random block matrix."""
    return cut_matrix(np.random.rand(sum(x_shapes), sum(y_shapes)), x_shapes, y_shapes)


def full_like(A, value, dtype=np.float64):
    """A matrix of the same kind and shape as A but filled with a single value."""
    if isinstance(A, BlockMatrix):
        new_matrix = []
        for i in range(A._stored_nb_blocks[0]):
            line = []
            for j in range(A._stored_nb_blocks[1]):
                line.append(full_like(A._stored_blocks[i, j], value, dtype=dtype))
            new_matrix.append(line)
        return A.__class__(new_matrix)
    elif isinstance(A, LowRankMatrix):
        return LowRankMatrix(np.ones((A.shape[0], 1)), np.full((1, A.shape[1]), value))
    elif isinstance(A, np.ndarray):
        return np.full_like(A, value, dtype=dtype)


def zeros_like(A, dtype=np.float64):
    """A matrix of the same kind and shape as A but filled with zeros."""
    return full_like(A, 0.0, dtype=dtype)


def ones_like(A, dtype=np.float64):
    """A matrix of the same kind and shape as A but filled with ones."""
    return full_like(A, 1.0, dtype=dtype)


def identity_like(A, dtype=np.float64):
    """A identity matrix of the same kind and shape as A."""
    if isinstance(A, BlockMatrix):
        I = []
        for i in range(A._stored_nb_blocks[0]):
            line = []
            for j in range(A._stored_nb_blocks[1]):
                if i == j:
                    line.append(identity_like(A._stored_blocks[i, j], dtype=dtype))
                else:
                    line.append(zeros_like(A._stored_blocks[i, j], dtype=dtype))
            I.append(line)
        return A.__class__(I)
    elif isinstance(A, np.ndarray):
        return np.eye(A.shape[0], A.shape[1], dtype=dtype)
