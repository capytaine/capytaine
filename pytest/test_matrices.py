#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.Toeplitz_matrices import *


def test_BlockToeplitz():
    A = BlockToeplitzMatrix([np.array([[i]]) for i in range(5)])
    assert A.nb_blocks == 5
    assert A.block_size == 1
    assert np.all(A == A.full_matrix())
    assert np.all(A == np.array([[0, 1, 2, 3, 4],
                                 [1, 0, 1, 2, 3],
                                 [2, 1, 0, 1, 2],
                                 [3, 2, 1, 0, 1],
                                 [4, 3, 2, 1, 0]])
                  )
    assert isinstance(A + 0.1, BlockToeplitzMatrix)
    assert isinstance(0.1 + A, BlockToeplitzMatrix)
    assert isinstance(A - 1, BlockToeplitzMatrix)
    assert isinstance(1 - A, BlockToeplitzMatrix)
    assert np.all(A - 1 == np.array([[-1, 0, 1, 2, 3],
                                     [0, -1, 0, 1, 2],
                                     [1, 0, -1, 0, 1],
                                     [2, 1, 0, -1, 0],
                                     [3, 2, 1, 0, -1]])
                  )
    assert isinstance(A * 2, BlockToeplitzMatrix)
    assert isinstance(3 * A, BlockToeplitzMatrix)
    assert isinstance(A / 4, BlockToeplitzMatrix)
    assert isinstance(5/(A+1), BlockToeplitzMatrix)

    B = np.random.rand(5, 5)
    assert np.all(A + B == A.full_matrix() + B)
    assert np.all(A * B == A.full_matrix() * B)
    assert np.all(A @ B == A.full_matrix() @ B)

    C = np.random.rand(30, 30)
    Z = BlockToeplitzMatrix([np.zeros((3, 3)) for _ in range(10)])
    assert Z.nb_blocks == 10
    assert Z.block_size == 3
    assert np.all(Z.full_matrix() == np.zeros((30, 30)))
    assert np.all(2 * Z == Z)
    assert np.all(Z * 2 == Z)
    assert np.all(Z * C == Z.full_matrix())
    # assert np.all(C * Z == Z.full_matrix())
    assert np.all(Z @ C == Z.full_matrix())
    # assert np.all(C @ Z == Z.full_matrix())

def test_solve():
    # 2x2 blocks
    A1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A2 = np.array([[5, 4, 2], [8, 0, 1], [6, 7, 3]])
    A = BlockToeplitzMatrix([A1, A2])
    b = np.array([3, 5, 7, 11, 13, 17])

    x_toe = solve(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_toe, x_dumb, rtol=1e-6)

