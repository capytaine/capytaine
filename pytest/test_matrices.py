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
    assert A.astype(np.complex128).full_matrix().dtype == np.complex128
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

    x = np.array([10, 11, 12, 13 ,14])
    assert np.all(A @ x == A.full_matrix() @ x)
    assert np.all(A @ A.full_matrix() == A.full_matrix() @ A.full_matrix())

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


def test_Circulant():
    matrix_list = [np.array([[i]]) for i in range(4)]
    A = BlockCirculantMatrix(matrix_list)
    assert A.nb_blocks == 6
    assert A.block_size == 1
    assert np.all(A == A.full_matrix())
    assert np.all(A == np.array([[0, 1, 2, 3, 2, 1],
                                 [1, 0, 1, 2, 3, 2],
                                 [2, 1, 0, 1, 2, 3],
                                 [3, 2, 1, 0, 1, 2],
                                 [2, 3, 2, 1, 0, 1],
                                 [1, 2, 3, 2, 1, 0],
                                 ])
                  )

    B = BlockCirculantMatrix(matrix_list, size=6)
    assert np.all(A == B)

    C = BlockCirculantMatrix(matrix_list, size=5)
    assert np.all(C == np.array([[0, 1, 2, 2, 1],
                                 [1, 0, 1, 2, 2],
                                 [2, 1, 0, 1, 2],
                                 [2, 2, 1, 0, 1],
                                 [1, 2, 2, 1, 0],
                                 ])
                  )

    assert isinstance(A + 0.1, BlockCirculantMatrix)
    assert (A + 0.1).nb_blocks == A.nb_blocks
    assert (A + 0.1).block_size == A.block_size
    assert isinstance(0.1 + A, BlockCirculantMatrix)
    assert isinstance(A - 1, BlockCirculantMatrix)
    assert isinstance(1 - A, BlockCirculantMatrix)
    assert isinstance(A * 2, BlockCirculantMatrix)
    assert (A * 2).nb_blocks == A.nb_blocks
    assert (A * 2).block_size == A.block_size
    assert isinstance(3 * A, BlockCirculantMatrix)
    assert isinstance(A / 4, BlockCirculantMatrix)
    assert (A / 4).nb_blocks == A.nb_blocks
    assert (A / 4).block_size == A.block_size
    assert isinstance(5/(A+1), BlockCirculantMatrix)

    I = block_circulant_identity(A.nb_blocks, A.block_size)
    assert isinstance(A + I, BlockCirculantMatrix)

    I2 = block_circulant_identity(C.nb_blocks, C.block_size)
    assert isinstance(C + I2, BlockCirculantMatrix)


def test_solve_2x2():
    # 2x2 blocks
    A1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A2 = np.array([[5, 4, 2], [8, 0, 1], [6, 7, 3]])
    A = BlockToeplitzMatrix([A1, A2])

    b = np.random.rand(6)

    x_toe = solve(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_toe, x_dumb, rtol=1e-6)


def test_solve_block_circulant():
    # Block Circulant Matrix
    A1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A2 = np.array([[5, 4, 2], [8, 0, 1], [6, 7, 3]])
    A3 = np.array([[0, 0, 3], [9, 3, 5], [7, 5, 6]])

    A = BlockCirculantMatrix([A1, A2, A3])
    b = np.random.rand(12)

    x_toe = solve(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_toe, x_dumb, rtol=1e-6)
