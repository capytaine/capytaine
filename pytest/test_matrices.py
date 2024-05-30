"""Test of the matrices submodule."""

import pytest

import numpy as np
from numpy.linalg import norm, matrix_rank

from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.block_toeplitz import *
from capytaine.matrices.builders import *
from capytaine.matrices.low_rank import LowRankMatrix, NoConvergenceOfACA

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

@pytest.fixture
def two_by_two_block_identity():
    A = BlockMatrix([
        [np.eye(2, 2), np.zeros((2, 2))],
        [np.zeros((2, 2)), np.eye(2, 2)]
    ])
    return A


def test_block_matrix_representation_of_identity_two_by_two(two_by_two_block_identity):
    # 2x2 block representation of the identity matrix
    A = two_by_two_block_identity
    assert A.shape == (4, 4)
    assert A.nb_blocks == (2, 2)
    assert A.block_shapes == ([2, 2], [2, 2])
    assert list(A._stored_block_positions()) == [[(0, 0)], [(0, 2)], [(2, 0)], [(2, 2)]]
    assert A.str_shape == "2×2×[2×2]"
    assert A.density == 1.0

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()
    b = np.random.rand(4)
    assert (A @ b == b).all()
    assert (A @ A == A).all()

    assert (A == identity_like(A)).all()
    assert (A.full_matrix() == np.eye(4, 4)).all()
    assert (np.array(A) == A.full_matrix()).all()


@pytest.mark.skipif(plt is None,
                     reason='matplotlib is not installed')
def test_block_matrix_representation_of_identity_patches(two_by_two_block_identity):
    A = two_by_two_block_identity
    patches = A._patches(global_frame=(10, 10))
    exp = {(10, 10), (12, 10), (10, 12), (12, 12)}
    assert {rectangle.get_xy() for rectangle in patches} == exp


@pytest.fixture
def block_three_rect(two_by_two_block_identity):
    A = two_by_two_block_identity
    B = BlockMatrix([[A, np.zeros((4, 1))], [np.zeros((1, 4)), np.ones((1, 1))]])
    return B


def test_block_matrix_representation_of_identity_rect(block_three_rect):
    # 2x2 block matrix containing one another block matrix and three regular matrices of different shapes
    B = block_three_rect
    assert B.shape == (5, 5)
    assert B.block_shapes == ([4, 1], [4, 1])
    assert (B.full_matrix() == np.eye(5, 5)).all()
    assert B.str_shape == "[[2×2×[2×2], 4×1], [1×4, 1×1]]"
    assert B.block_shapes == ([4, 1], [4, 1])
    assert list(B._stored_block_positions()) == [[(0, 0)], [(0, 4)], [(4, 0)], [(4, 4)]]


@pytest.mark.skipif(plt is None,
                     reason='matplotlib is not installed')
def test_block_matrix_representation_of_identity_rect_patch(block_three_rect):
    B = block_three_rect
    patches = B._patches(global_frame=(10, 10))
    exp = {(10, 10), (12, 10), (10, 12), (12, 12), (14, 10), (10, 14), (14, 14)}
    assert {rectangle.get_xy() for rectangle in patches} == exp

def test_block_diff_shapes():
    # 3x3 block matrix with blocks of different shapes
    C = random_block_matrix([1, 2, 4], [1, 2, 4])
    assert C.nb_blocks == (3, 3)
    assert C.block_shapes == ([1, 2, 4], [1, 2, 4])
    assert (ones_like(C).full_matrix() ==  np.ones(C.shape)).all()

    assert (cut_matrix(C.full_matrix(), *C.block_shapes) == C).all()

    assert (C @ random_block_matrix(C.block_shapes[1], [6])).block_shapes == ([1, 2, 4], [6])

    b = np.random.rand(7)
    assert np.allclose(C @ b, C.full_matrix() @ b)



def test_sparse_storage_of_block_toeplitz_matrices():
    A = BlockToeplitzMatrix(
        [[np.array([[i]]) for i in range(5)]]
    )
    assert np.all(A.full_matrix() ==
                  np.array([[0, 1, 2],
                            [4, 0, 1],
                            [3, 4, 0]]))
    assert A.density == 5 / 9

    assert not isinstance(A.no_toeplitz(), BlockToeplitzMatrix)
    assert np.all(A.no_toeplitz().full_matrix() == A.full_matrix())

    B = BlockSymmetricToeplitzMatrix(
        [[np.array([[i]]) for i in range(4)]]
    )
    assert np.all(B.full_matrix() ==
                  np.array([[0, 1, 2, 3],
                            [1, 0, 1, 2],
                            [2, 1, 0, 1],
                            [3, 2, 1, 0],
                            ]))
    assert B.density == 4 / 16

    C = BlockCirculantMatrix(
        [[np.array([[i]]) for i in range(4)]]
    )
    assert np.all(C.full_matrix() ==
                  np.array([[0, 1, 2, 3],
                            [3, 0, 1, 2],
                            [2, 3, 0, 1],
                            [1, 2, 3, 0],
                            ]))
    assert C.density == 4 / 16

    D = EvenBlockSymmetricCirculantMatrix(
        [[np.array([[i]]) for i in range(3)]]
    )
    assert np.all(D.full_matrix() ==
                  np.array([[0, 1, 2, 1],
                            [1, 0, 1, 2],
                            [2, 1, 0, 1],
                            [1, 2, 1, 0],
                            ]))
    assert D.density == 3 / 16

    E = OddBlockSymmetricCirculantMatrix(
        [[np.array([[i]]) for i in range(3)]]
    )
    assert np.all(E.full_matrix() ==
                  np.array([[0, 1, 2, 2, 1],
                            [1, 0, 1, 2, 2],
                            [2, 1, 0, 1, 2],
                            [2, 2, 1, 0, 1],
                            [1, 2, 2, 1, 0],
                            ]))
    assert E.density == 3 / 25


def test_block_toeplitz_matrices():
    # 2x2 block Toeplitz representation of the identity matrix
    A = BlockToeplitzMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.block_shapes == ([2, 2], [2, 2])
    assert A.block_shape == (2, 2)
    assert A.shape == (4, 4)
    assert list(A._stored_block_positions()) == [[(0, 0), (2, 2)], [(0, 2)], [(2, 0)]]
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    # Removing the Toeplitz structure.
    Ant = A.no_toeplitz()
    assert not isinstance(Ant, BlockToeplitzMatrix)
    assert Ant.all_blocks[0, 0] is Ant.all_blocks[1, 1]
    # When using no_toeplitz, the resulting BlockMatrix still contains several references to the same array object.

    from copy import deepcopy
    Ant_cp = deepcopy(Ant)
    assert Ant_cp.all_blocks[0, 0] is not Ant_cp.all_blocks[1, 1]
    # The deepcopy has made different copies of the blocks in the block matrix without explicit Toeplitz structure.

    # However,
    A_cp = deepcopy(A)
    assert A_cp.all_blocks[0, 0] is A_cp.all_blocks[1, 1]
    # The matrix with explicit Toeplitz structure still keep the structure when deepcopied.

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    assert isinstance(A.circulant_super_matrix, BlockCirculantMatrix)
    assert A.circulant_super_matrix.shape == (6, 6)

    b = np.random.rand(A.shape[0])
    assert np.allclose(A @ b, b)
    assert np.allclose(A.rmatvec(b), b)

    # 3x3 block random Toeplitz matrix
    B = BlockToeplitzMatrix([
        [
            random_block_matrix([1], [1]),
            random_block_matrix([1], [1]),
            random_block_matrix([1], [1]),
            random_block_matrix([1], [1]),
            random_block_matrix([1], [1]),
        ]
    ])
    assert B.shape == (3, 3)
    assert B.nb_blocks == (3, 3)
    assert B.block_shape == (1, 1)
    assert list(B._stored_block_positions()) == [[(0, 0), (1, 1), (2, 2)], [(0, 1), (1, 2)], [(0, 2)], [(2, 0)], [(1, 0), (2, 1)]]

    assert isinstance(B.circulant_super_matrix, BlockCirculantMatrix)
    assert B.circulant_super_matrix.shape == (5, 5)

    b = np.random.rand(B.shape[0])
    assert np.allclose(BlockMatrix.matvec(B, b), B.full_matrix() @ b)
    assert np.allclose(B @ b, B.full_matrix() @ b)
    assert np.allclose(B.rmatvec(b), b @ B.full_matrix())


def test_block_symmetric_toeplitz_matrices():
    # 2x2 block symmetric Toeplitz representation of the identity matrix
    A = BlockSymmetricToeplitzMatrix([
        [np.eye(2, 2), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (2, 2)
    assert A.block_shapes == ([2, 2], [2, 2])
    assert A.shape == (4, 4)
    assert A.block_shape == (2, 2)
    assert list(A._stored_block_positions()) == [[(0, 0), (2, 2)], [(0, 2), (2, 0)]]
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    b = np.random.rand(4)
    assert np.allclose(A @ b, b)
    assert np.allclose(A.rmatvec(b), b)

    # The same as a simpler block matrix
    A2 = BlockMatrix(A.all_blocks)
    assert (A2.full_matrix() == A.full_matrix()).all()

    # 2x2 blocks symmetric Toeplitz matrix of non-squared matrix
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

    # 2x2 block symmetrix Toeplitz matrix composed of block matrices
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

    list_of_matrices = [B, B + ones_like(B), B, B - ones_like(B)]
    BB_fft = BlockSymmetricToeplitzMatrix.fft_of_list(*list_of_matrices)
    full_BB_fft = np.fft.fft(np.array([A.full_matrix() for A in list_of_matrices]), axis=0)
    assert np.allclose(full_BB_fft, np.array([A.full_matrix() for A in BB_fft]))

    # 2x2 block symmetric Toeplitz matrix composed of the previous block symmetric Toeplitz matrices
    C = BlockSymmetricToeplitzMatrix([
        [A, B]
    ])
    assert C.shape == (8, 8)

    b = np.random.rand(8)
    assert np.allclose(C @ b, C.full_matrix() @ b)
    assert np.allclose(C.rmatvec(b), b @ C.full_matrix())

    # We need to go deeper
    D = BlockMatrix([
        [C, np.zeros(C.shape)]
    ])
    assert D.shape == (8, 16)

    b = np.random.rand(16)
    assert np.allclose(D @ b, D.full_matrix() @ b)


def test_block_circulant_matrix():
    # 5x5 block symmetric circulant matrix representation of the identity matrix
    A = BlockCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (3, 3)
    assert A.shape == (6, 6)
    assert A.block_shapes == ([2, 2, 2], [2, 2, 2])
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    b = np.random.rand(A.shape[0])
    assert np.allclose(A @ b, A.full_matrix() @ b)
    assert np.allclose(A.rmatvec(b), b @ A.full_matrix())

    # Nested matrix
    B = BlockCirculantMatrix([[A, 2*A, 3*A]])
    assert B.nb_blocks == (3, 3)
    assert B.shape == (18, 18)
    assert (B.all_blocks[0, 0] == A).all()

    C = BlockCirculantMatrix([[np.array([[i]]) for i in range(5)]])
    b = np.random.rand(C.shape[0])
    assert np.allclose(BlockMatrix.matvec(C, b), C.full_matrix() @ b)
    assert np.allclose(C @ b, C.full_matrix() @ b)
    assert np.allclose(C.rmatvec(b), b @ C.full_matrix())


def test_even_block_symmetric_circulant_matrix():
    # 5x5 block symmetric circulant matrix representation of the identity matrix
    A = EvenBlockSymmetricCirculantMatrix([
        [np.eye(2, 2), np.zeros((2, 2)), np.zeros((2, 2))]
    ])
    assert A.nb_blocks == (4, 4)
    assert A.shape == (8, 8)
    assert A.block_shapes == ([2, 2, 2, 2], [2, 2, 2, 2])
    assert (A.full_matrix() == np.eye(*A.shape)).all()

    assert ((A + A)/2 == A).all()
    assert (-A).min() == -1
    assert (2*A).max() == 2
    assert (A*A == A).all()

    b = np.random.rand(A.shape[0])
    assert np.allclose(A @ b, A.full_matrix() @ b)
    assert np.allclose(A.rmatvec(b), b @ A.full_matrix())

    B = EvenBlockSymmetricCirculantMatrix([
        [A, A, A]
    ])
    assert B.nb_blocks == (4, 4)
    assert B.shape == (32, 32)
    assert B.all_blocks[0, 0] is A
    assert (B.all_blocks[0, 1] == B.all_blocks[2, 3]).all()

    b = np.random.rand(B.shape[0])
    assert np.allclose(B @ b, B.full_matrix() @ b)
    assert np.allclose(B.rmatvec(b), b @ B.full_matrix())


def test_complex_valued_matrices():
    R = random_block_matrix([2, 2], [2, 2])
    I = random_block_matrix([2, 2], [2, 2])
    A = R + 1j * I
    assert np.allclose(A.full_matrix(), R.full_matrix() + 1j * I.full_matrix())

    B = EvenBlockSymmetricCirculantMatrix([[A.full_matrix(), 2*A.full_matrix()]])
    c = np.random.rand(B.shape[1]) + 1j * np.random.rand(B.shape[1])
    assert np.allclose(B @ c, B.full_matrix() @ c)


def test_low_rank_blocks():
    n = 10

    # Test initialization
    a, b = np.random.rand(n, 1), np.random.rand(1, n)
    LR = LowRankMatrix(a, b)
    assert LR.shape == LR.full_matrix().shape
    assert np.all(np.array(LR) == LR.full_matrix())
    assert matrix_rank(LR.full_matrix()) == LR.rank == 1
    assert LR.density == 2 * n / n ** 2

    a, b = np.random.rand(n, 2), np.random.rand(2, n)
    LR = LowRankMatrix(a, b)
    assert LR.shape == LR.full_matrix().shape
    assert matrix_rank(LR.full_matrix()) == LR.rank == 2

    # Test creation from SVD
    A = np.arange(n**2).reshape((n, n)) + 1.0
    dumb_low_rank = LowRankMatrix.from_full_matrix_with_SVD(A, n)
    assert np.allclose(dumb_low_rank.full_matrix() - A, 0.0)

    A_rank_1 = LowRankMatrix.from_full_matrix_with_SVD(A, 1)
    assert matrix_rank(A_rank_1.full_matrix()) == A_rank_1.rank == 1

    # Test recompression
    recompressed = dumb_low_rank.recompress(new_rank=2)
    assert recompressed.rank == matrix_rank(recompressed.full_matrix()) == 2

    recompressed = dumb_low_rank.recompress(tol=1e-1)
    assert recompressed.rank <= dumb_low_rank.rank

    # Test multiplication with vector
    b = np.random.rand(n)
    assert np.allclose(A_rank_1 @ b, A_rank_1.full_matrix() @ b)

    # Test creation with ACA
    full_A_rank_1 = A_rank_1.full_matrix()
    A_rank_1_again = LowRankMatrix.from_full_matrix_with_ACA(full_A_rank_1, max_rank=5)
    assert matrix_rank(A_rank_1_again.full_matrix()) == 1
    assert np.allclose(A_rank_1_again.full_matrix(), full_A_rank_1)

    # Test creation from function with ACA
    X = np.linspace(0, 1, n)
    Y = np.linspace(10, 11, n)

    def f(i, j):
        return 1/abs(X[i] - Y[j])

    S = np.array([[f(i, j) for j in range(n)] for i in range(n)])
    SLR = LowRankMatrix.from_function_with_ACA(f, n, n, max_rank=2, tol=1e-2)
    assert SLR.shape == (n, n)
    assert np.allclose(SLR.full_matrix(), S, atol=1e-2)

    with pytest.raises(NoConvergenceOfACA):
        LowRankMatrix.from_function_with_ACA(f, n, n, max_rank=1, tol=1e-3)

    summed = SLR + A_rank_1
    assert summed.rank == 1


def test_2_in_1_ACA_with_identical_matrices():
    n = 5
    A = np.fromfunction(lambda i, j: 1/abs(i - (100 + j)), (n, n))
    B = A.copy()

    def get_row_func(i):
        return A[i, :], B[i, :]

    def get_col_func(j):
        return A[:, j], B[:, j]

    lrA, lrB = LowRankMatrix.from_rows_and_cols_functions_with_multi_ACA(
        get_row_func, get_col_func, n, n, nb_matrices=2, tol=1e-3
    )

    assert (lrA.full_matrix() == lrB.full_matrix()).all()
    assert norm(lrA.full_matrix() - A, 'fro')/norm(A, 'fro') < 3e-2


def test_2_in_1_ACA_with_different_matrices():
    n = 5
    A = np.arange(1, 1+n**2).reshape((n, n)) + np.random.rand(n, n)
    B = A.T

    def get_row_func(i):
        return A[i, :], B[i, :]

    def get_col_func(j):
        return A[:, j], B[:, j]

    lrA, lrB = LowRankMatrix.from_rows_and_cols_functions_with_multi_ACA(
        get_row_func, get_col_func, n, n, nb_matrices=2, max_rank=3
    )

    assert matrix_rank(lrA.full_matrix()) == matrix_rank(lrB.full_matrix()) == 3


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

def test_rmatvec_lowrank():

    left_matrix = np.arange(6).reshape((3,2))
    right_matrix = np.arange(8).reshape((2,4))

    v = np.ones(3)

    M_lowrank = LowRankMatrix(left_matrix, right_matrix)
    M_full = left_matrix@right_matrix

    vM_lowrank = M_lowrank.__rmatmul__(v)
    vM_full = v @ M_full

    assert np.all(vM_lowrank==vM_full)

def test_access_block_by_path():
    A = BlockMatrix([[random_block_matrix([2, 2], [2, 2]), random_block_matrix([2, 2], [2, 2])],
                     [random_block_matrix([2, 2], [2, 2]), random_block_matrix([2, 2], [2, 2])]])
    assert A.access_block_by_path([0]) is A._stored_blocks[0, 0]
    assert A.access_block_by_path([0, 0]) is A._stored_blocks[0, 0]._stored_blocks[0, 0]
    assert A.access_block_by_path([0, 1]) is A._stored_blocks[0, 0]._stored_blocks[1, 1]
    assert A.access_block_by_path([1, 0]) is A._stored_blocks[1, 1]._stored_blocks[0, 0]
