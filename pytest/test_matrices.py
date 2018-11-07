#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.matrices.block_matrices import *
from capytaine.matrices.block_toeplitz_matrices import *
from capytaine.matrices.builders import *
from capytaine.matrices.low_rank_blocks import LowRankMatrix
from capytaine.matrices.linear_solvers import solve_directly


def test_block_matrices():
    A = BlockMatrix([
        [np.eye(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.eye(2, 2)]
    ])
    assert A.shape == (4, 4)
    assert A.nb_blocks == (2, 2)
    assert A.block_shapes == ([2, 2], [2, 2])
    assert list(A._stored_block_positions()) == [[(0, 0)], [(0, 2)], [(2, 0)], [(2, 2)]]
    assert repr(A) == "BlockMatrix(nb_blocks=(2, 2), shape=(4, 4))"

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()
    b = np.random.rand(4)
    assert (A @ b == b).all()
    assert (A @ A == A).all()

    patches = A._patches(global_frame=(10, 10))
    assert {rectangle.get_xy() for rectangle in patches} == {(10, 10), (12, 10), (10, 12), (12, 12)}

    assert (A.T == A).all()
    assert not (A.T != A).any()
    assert (A == identity_like(A)).all()
    assert (A.full_matrix() == np.eye(4, 4)).all()

    B = BlockMatrix([[A, np.zeros((4, 1))], [np.zeros((1, 4)), np.ones((1, 1))]])
    assert B.shape == (5, 5)
    assert B.block_shapes == ([4, 1], [4, 1])
    assert (B == B.T).all()
    assert (B.full_matrix() == np.eye(5, 5)).all()
    assert repr(B) == "BlockMatrix(nb_blocks=(2, 2), shape=(5, 5))"
    assert B.block_shapes == ([4, 1], [4, 1])
    assert list(B._stored_block_positions()) == [[(0, 0)], [(0, 4)], [(4, 0)], [(4, 4)]]

    patches = B._patches(global_frame=(10, 10))
    assert {rectangle.get_xy() for rectangle in patches} == {(10, 10), (12, 10), (10, 12), (12, 12),
                                                             (14, 10), (10, 14), (14, 14)}

    C = random_block_matrix([1, 2, 4], [1, 2, 4])
    assert C.nb_blocks == (3, 3)
    assert C.block_shapes == ([1, 2, 4], [1, 2, 4])
    assert (ones_like(C).full_matrix() ==  np.ones(C.shape)).all()

    assert (cut_matrix(C.full_matrix(), *C.block_shapes) == C).all()

    assert (C @ random_block_matrix(C.block_shapes[1], [6])).block_shapes == ([1, 2, 4], [6])

    b = np.random.rand(7)
    assert np.allclose(C @ b, C.full_matrix() @ b)


def test_block_toeplitz_matrices():
    A = BlockSymmetricToeplitzMatrix([
        [np.eye(2, 2), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (2, 2)
    assert A.block_shapes == ([2, 2], [2, 2])
    assert A.shape == (4, 4)
    assert A.block_shape == (2, 2)
    assert list(A._stored_block_positions()) == [[(0, 0), (2, 2)], [(0, 2), (2, 0)]]
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert (A._index_grid() == np.array([[0, 1], [1, 0]])).all()
    assert (A == A.T).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    b = np.random.rand(4)
    assert np.allclose(A @ b, b)

    A2 = BlockMatrix(A.all_blocks)
    assert (A2.full_matrix() == A.full_matrix()).all()

    A3 = BlockSymmetricToeplitzMatrix([
        [np.ones((3, 1)), np.zeros((3, 1))]
    ])
    assert A3.nb_blocks == (2, 2)
    assert A3.block_shapes == ([3, 3], [1, 1])
    assert A3.shape == (6, 2)
    assert A3.block_shape == (3, 1)

    with pytest.raises(AssertionError):
        BlockSymmetricToeplitzMatrix([
            [np.ones((1, 1)), np.zeros((2, 2))]
        ])

    B = BlockSymmetricToeplitzMatrix([
        [random_block_matrix([1, 1], [1, 1]), random_block_matrix([1, 1], [1, 1])]
    ])
    assert B.nb_blocks == (2, 2)
    assert B._stored_nb_blocks == (1, 2)
    assert B.block_shapes == ([2, 2], [2, 2])
    assert B.block_shape == (2, 2)
    assert B.shape == (4, 4)

    assert isinstance(B.all_blocks[0, 0], BlockMatrix)
    assert (B.all_blocks[0, 0] == B.all_blocks[1, 1]).all()

    C = BlockSymmetricToeplitzMatrix([
        [A, B]
    ])
    assert C.shape == (8, 8)

    b = np.random.rand(8)
    assert np.allclose(C @ b, C.full_matrix() @ b)

    D = BlockMatrix([
        [C, np.zeros(C.shape)]
    ])
    assert D.shape == (8, 16)

    b = np.random.rand(16)
    assert np.allclose(D @ b, D.full_matrix() @ b)


def test_even_block_circulant_matrix():
    A = EvenBlockSymmetricCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (4, 4)
    assert A.shape == (8, 8)
    assert A.block_shapes == ([2, 2, 2, 2], [2, 2, 2, 2])
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

    b = np.random.rand(8)
    assert np.allclose(A @ b, A.full_matrix() @ b)

    B = EvenBlockSymmetricCirculantMatrix([
        [A, A, A]
    ])
    assert B.nb_blocks == (4, 4)
    assert B.shape == (32, 32)
    assert B.all_blocks[0, 0] is A
    assert (B.all_blocks[0, 1] == B.all_blocks[2, 3]).all()

    b = np.random.rand(32)
    assert np.allclose(B @ b, B.full_matrix() @ b)


def test_odd_block_circulant_matrix():
    A = OddBlockSymmetricCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (5, 5)
    assert A.shape == (10, 10)
    assert A.block_shapes == ([2, 2, 2, 2, 2], [2, 2, 2, 2, 2])
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert (A._index_grid() == np.array([
        [0, 1, 2, 2, 1], [1, 0, 1, 2, 2], [2, 1, 0, 1, 2], [2, 2, 1, 0, 1], [1, 2, 2, 1, 0]
    ])).all()
    assert (A == A.T).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    assert A.first_block_line.ndim == 3
    assert A.first_block_line.shape[:1] == (5,)

    b = np.random.rand(10)
    assert np.allclose(A @ b, A.full_matrix() @ b)

    B = OddBlockSymmetricCirculantMatrix([
        [A, A, A]
    ])
    assert B.nb_blocks == (5, 5)
    assert B.shape == (50, 50)
    assert B.all_blocks[0, 0] is A
    assert (B.all_blocks[0, 1] == B.all_blocks[2, 3]).all()

    b = np.random.rand(50)
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
    A = EvenBlockSymmetricCirculantMatrix([
        [np.random.rand(3, 3) for _ in range(4)]
    ])
    b = np.random.rand(A.shape[0])

    x_circ = solve_directly(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_circ, x_dumb, rtol=1e-6)

    A = EvenBlockSymmetricCirculantMatrix([
        [random_block_matrix([1, 1], [1, 1]) for _ in range(5)]
    ])
    b = np.random.rand(A.shape[0])

    x_circ = solve_directly(A, b)
    x_dumb = np.linalg.solve(A.full_matrix(), b)

    assert np.allclose(x_circ, x_dumb, rtol=1e-6)


def test_low_rank_blocks():
    n = 10

    # Test initialization
    a, b = np.random.rand(n, 1), np.random.rand(1, n)
    LR = LowRankMatrix(a, b)
    assert LR.shape == LR.full_matrix().shape
    assert np.linalg.matrix_rank(LR.full_matrix()) == LR.rank == 1

    a, b = np.random.rand(n, 2), np.random.rand(2, n)
    LR = LowRankMatrix(a, b)
    assert LR.shape == LR.full_matrix().shape
    assert np.linalg.matrix_rank(LR.full_matrix()) == LR.rank == 2

    # Test creation from SVD
    A = np.arange(n**2).reshape((n, n)) + np.random.rand(n, n)
    dumb_low_rank = LowRankMatrix.from_full_matrix_with_SVD(A, n)
    assert np.allclose(dumb_low_rank.full_matrix() - A, 0.0)

    A_rank_1 = LowRankMatrix.from_full_matrix_with_SVD(A, 1)
    assert np.linalg.matrix_rank(A_rank_1.full_matrix()) == A_rank_1.rank == 1

    # Test recompression
    recompressed = dumb_low_rank.recompress(new_rank=2)
    assert recompressed.rank == np.linalg.matrix_rank(recompressed.full_matrix()) == 2

    recompressed = dumb_low_rank.recompress(tol=1e-1)
    assert recompressed.rank <= dumb_low_rank.rank

    # Test multiplication with vector
    b = np.random.rand(n)
    assert np.allclose(A_rank_1 @ b, A_rank_1.full_matrix() @ b)

    # Test creation with ACA
    full_A_rank_1 = A_rank_1.full_matrix()
    A_rank_1_again = LowRankMatrix.from_full_matrix_with_ACA(full_A_rank_1, max_rank=5)
    assert np.linalg.matrix_rank(A_rank_1_again.full_matrix()) == 1
    assert np.allclose(A_rank_1_again.full_matrix(), full_A_rank_1)

    # Test creation from function with ACA
    X = np.linspace(0, 1, n)
    Y = np.linspace(5, 6, n)

    def f(i, j):
        return 1/abs(X[i] - Y[j])

    S = np.array([[f(i, j) for j in range(n)] for i in range(n)])
    SLR = LowRankMatrix.from_function_with_ACA(f, n, n, max_rank=2, tol=1e-5)
    assert SLR.shape == (n, n)
    assert np.allclose(SLR.full_matrix(), S, atol=1e-4)

    summed = SLR + A_rank_1
    assert summed.rank == 1


def test_hierarchical_matrix():
    n = 30
    X = np.linspace(0, 1, n)
    Y = np.linspace(10, 11, n)

    def f(i, j):
        return 1/abs(X[i] - Y[j])

    S = np.array([[f(i, j) for j in range(n)] for i in range(n)])

    HS = BlockMatrix(
        [
            [S[:n//2, :n//2], LowRankMatrix.from_full_matrix_with_ACA(S[n//2:, :n//2], tol=1e-5)],
            [LowRankMatrix.from_full_matrix_with_ACA(S[:n//2, n//2:], tol=1e-5), S[n//2:, n//2:]],
         ]
    )

    assert np.allclose(HS.full_matrix(), S, rtol=2e-1)

    doubled = HS + HS
    assert np.allclose(2*S, doubled.full_matrix(), rtol=2e-1)


