#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.Toeplitz_matrices import *


def test_BlockToeplitz():
    A = BlockToeplitzMatrix([np.array([[i]]) for i in range(5)])
    assert A.nb_blocks == 5
    assert A.block_size == 1
    assert np.all(A.full_matrix() == np.array([[0, 1, 2, 3, 4],
                                               [1, 0, 1, 2, 3],
                                               [2, 1, 0, 1, 2],
                                               [3, 2, 1, 0, 1],
                                               [4, 3, 2, 1, 0]])
                  )

    A = BlockToeplitzMatrix([np.zeros((3, 3)) for _ in range(10)])
    assert A.nb_blocks == 10
    assert A.block_size == 3
    assert np.all(A.full_matrix() == np.zeros((30, 30)))

def test_solve():
    # 2x2 blocks
    A1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A2 = np.array([[5, 4, 2], [8, 0, 1], [6, 7, 3]])
    A = BlockToeplitzMatrix([A1, A2])
    b = np.array([3, 5, 7, 11, 13, 17])

    x_toe = solve(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_toe, x_dumb, rtol=1e-6)

