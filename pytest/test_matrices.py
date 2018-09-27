#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.matrices.block_matrices import *
from capytaine.matrices.block_toeplitz_matrices import *
from capytaine.matrices.builders import *
from capytaine.matrices.linear_solvers import solve_directly


def test_block_matrices():
    A = BlockMatrix([
        [np.eye(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.eye(2, 2)]
    ])
    assert A.shape == (4, 4)
    assert A.nb_blocks == (2, 2)
    assert set(A._block_positions_list) == {(0, 0), (2, 0), (0, 2), (2, 2)}
    assert repr(A) == "BlockMatrix(nb_blocks=(2, 2), shape=(4, 4))"

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()
    b = np.random.rand(4)
    assert (A @ b == b).all()
    assert (A @ A == A).all()

    patches = A._patches(global_shift=(10, 10))
    assert {rectangle.get_xy() for rectangle in patches} == {(10, 10), (12, 10), (10, 12), (12, 12)}

    assert (A.T == A).all()
    assert not (A.T != A).any()
    assert (A == identity_like(A)).all()
    assert (A.full_matrix() == np.eye(4, 4)).all()

    B = BlockMatrix([[A, np.zeros((4, 1))], [np.zeros((1, 4)), np.ones((1, 1))]])
    assert B.shape == (5, 5)
    assert (B == B.T).all()
    assert (B.full_matrix() == np.eye(5, 5)).all()
    assert repr(B) == "BlockMatrix(nb_blocks=(2, 2), shape=(5, 5))"
    assert B.block_shapes == ([4, 1], [4, 1])

    patches = B._patches(global_shift=(10, 10))
    assert {rectangle.get_xy() for rectangle in patches} == {(10, 10), (12, 10), (10, 12), (12, 12),
                                                             (14, 10), (10, 14), (14, 14)}

    C = random_block_matrix([1, 2, 4], [1, 2, 4])
    assert C.nb_blocks == (3, 3)
    assert C.block_shapes == ([1, 2, 4], [1, 2, 4])
    assert (ones_like(C).full_matrix() ==  np.ones(C.shape)).all()

    assert (cut_matrix(C.full_matrix(), *C.block_shapes) == C).all()

    assert (C @ random_block_matrix(C.block_shapes[1], [6])).block_shapes == ([1, 2, 4], [6])


def test_block_toeplitz_matrices():
    A = BlockSymmetricToeplitzMatrix([
        [np.eye(2, 2), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (2, 2)
    assert A.shape == (4, 4)
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert (A._index_grid() == np.array([[0, 1], [1, 0]])).all()
    assert (A == A.T).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    b = np.random.rand(4)
    assert (A @ b == b).all()

    A2 = BlockMatrix(A.all_blocks)
    assert (A2.full_matrix() == A.full_matrix()).all()

    B = BlockSymmetricToeplitzMatrix([
        [random_block_matrix([1, 1], [1, 1]), random_block_matrix([1, 1], [1, 1])]
    ])
    assert B.nb_blocks == (2, 2)
    assert B._nb_stored_blocks == (1, 2)
    assert B.block_shapes == ([2, 2], [2, 2])
    assert B.block_shape == (2, 2)
    assert B.shape == (4, 4)

    assert isinstance(B.all_blocks[0, 0], BlockMatrix)
    assert (B.all_blocks[0, 0] == B.all_blocks[1, 1]).all()

    C = BlockSymmetricToeplitzMatrix([
        [A, B]
    ])
    assert C.shape == (8, 8)

    D = BlockMatrix([
        [C, np.zeros(C.shape)]
    ])
    assert D.shape == (8, 16)

    b = np.random.rand(16)
    assert np.allclose(D @ b, D.full_matrix() @ b)


def test_block_circulant_matrix():
    A = BlockSymmetricCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (4, 4)
    assert A.shape == (8, 8)
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert (A._index_grid() == np.array([
        [0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]
    ])).all()
    assert (A == A.T).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    assert A.first_block_line.ndim == 3
    assert A.first_block_line.shape[:1] == (4,)

    A2 = BlockSymmetricCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ], size=5)
    assert (A2.full_matrix() == np.eye(*A2.shape)).all()

    B = BlockSymmetricCirculantMatrix([
        [A, A, A]
    ])
    assert B.nb_blocks == (4, 4)
    assert B.shape == (32, 32)
    assert B.all_blocks[0, 0] is A
    assert (B.all_blocks[0, 1] == B.all_blocks[2, 3]).all()

    b = np.random.rand(32)
    assert np.allclose(B @ b, B.full_matrix() @ b)


def test_solve_2x2():
    # 2x2 blocks
    A = BlockSymmetricToeplitzMatrix([
        [np.random.rand(3, 3) for _ in range(2)]
    ])
    b = np.random.rand(A.shape[0])

    x_toe = solve_directly(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_toe, x_dumb, rtol=1e-6)


def test_solve_block_circulant():
    # Block Circulant Matrix
    A = BlockSymmetricCirculantMatrix([
        [np.random.rand(3, 3) for _ in range(4)]
    ])
    b = np.random.rand(A.shape[0])

    x_circ = solve_directly(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_circ, x_dumb, rtol=1e-6)

    A = BlockSymmetricCirculantMatrix([
        [random_block_matrix([1, 1], [1, 1]) for _ in range(5)]
    ])
    b = np.random.rand(A.shape[0])

    x_circ = solve_directly(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_circ, x_dumb, rtol=1e-6)
