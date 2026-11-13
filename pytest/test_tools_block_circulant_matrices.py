# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import gc
import weakref

import numpy as np
import pytest
from capytaine.tools.block_circulant_matrices import (
        BlockCirculantMatrix, NestedBlockCirculantMatrix, lu_decompose,
        )

RNG = np.random.default_rng(seed=0)

def _rand(*size):
    return RNG.normal(size=size) + 1j*RNG.normal(size=size)


def test_2x2_block_circulant_matrices():
    A = BlockCirculantMatrix([
        2*np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
    ])
    full_A = np.array(A)
    assert full_A.shape == (4, 4)

    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
            A @ b,
            full_A @ b,
            )
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )


def test_deeper_2x2_block_circulant_matrices():
    A = BlockCirculantMatrix([
        _rand(2, 2, 3),
        _rand(2, 2, 3),
    ])
    full_A = np.array(A)
    assert full_A.shape == (4, 4, 3)


def test_3x3_block_circulant_matrices():
    A = BlockCirculantMatrix([
        np.eye(2) + _rand(2, 2),
        2*np.eye(2) + _rand(2, 2),
        3*np.eye(2) + _rand(2, 2)
    ])
    full_A = np.array(A)
    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
            A @ b,
            full_A @ b,
            )
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )

def test_4x4_block_circulant_matrices():
    A = BlockCirculantMatrix([
        np.eye(2) + _rand(2, 2),
        2*np.eye(2) + _rand(2, 2),
        3*np.eye(2) + _rand(2, 2),
        4*np.eye(2) + _rand(2, 2)
    ])
    full_A = np.array(A)
    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
            A @ b,
            full_A @ b,
            )
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )

def test_10x10_block_circulant_matrices():
    A = BlockCirculantMatrix([
        (lambda: np.eye(2) + _rand(2, 2))()
        for _ in range(10)
    ])
    full_A = np.array(A)
    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
            A @ b,
            full_A @ b,
            )
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )


def test_2x2x2x2_nested_block_circulant_matrices():
    """Test NestedBlockCirculantMatrix with 4 blocks (2x2 outer structure)"""
    A = NestedBlockCirculantMatrix([
        2*np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
    ])
    full_A = np.array(A)
    assert full_A.shape == (8, 8)

    b = RNG.normal(size=(A.shape[0],)) + 1j*RNG.normal(size=(A.shape[0],))
    assert np.allclose(
        A @ b,
        full_A @ b,
    )
    assert np.allclose(
        A.solve(b),
        np.linalg.solve(full_A, b)
    )

    # Test that the conversion to BlockCirculantMatrix produces equivalent matrix
    bc = A._block_circulant_matrix
    assert isinstance(bc, BlockCirculantMatrix)
    assert np.allclose(np.array(bc), full_A)

    # Test block_diagonalize
    bd = A.block_diagonalize()
    assert bd.nb_blocks == 4  # 2x2x2x2 fully diagonalized -> 4 blocks


def test_3x2_nested_block_circulant_matrices():
    """Test NestedBlockCirculantMatrix with 6 blocks (3x2 structure)"""
    A = NestedBlockCirculantMatrix([
        2*np.eye(3) + _rand(3, 3),
        np.eye(3) + _rand(3, 3),
        np.eye(3) + _rand(3, 3),
        np.eye(3) + _rand(3, 3),
        np.eye(3) + _rand(3, 3),
        np.eye(3) + _rand(3, 3),
    ])
    full_A = np.array(A)
    assert full_A.shape == (18, 18)

    b = RNG.normal(size=(A.shape[0],)) + 1j*RNG.normal(size=(A.shape[0],))
    assert np.allclose(
        A @ b,
        full_A @ b,
    )
    assert np.allclose(
        A.solve(b),
        np.linalg.solve(full_A, b)
    )

    # Test that the conversion to BlockCirculantMatrix produces equivalent matrix
    bc = A._block_circulant_matrix
    assert isinstance(bc, BlockCirculantMatrix)
    assert np.allclose(np.array(bc), full_A)


def test_4x2_nested_block_circulant_matrices():
    """Test NestedBlockCirculantMatrix with 8 blocks (4x2 structure)"""
    A = NestedBlockCirculantMatrix([
        2*np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
    ])
    full_A = np.array(A)
    assert full_A.shape == (16, 16)

    b = RNG.normal(size=(A.shape[0],)) + 1j*RNG.normal(size=(A.shape[0],))
    assert np.allclose(
        A @ b,
        full_A @ b,
    )
    assert np.allclose(
        A.solve(b),
        np.linalg.solve(full_A, b)
    )

    # Test that the conversion to BlockCirculantMatrix produces equivalent matrix
    bc = A._block_circulant_matrix
    assert isinstance(bc, BlockCirculantMatrix)
    assert np.allclose(np.array(bc), full_A)


def test_nested_block_circulant_matrix_structure():
    """Test that NestedBlockCirculantMatrix matches expected structure from docstring"""
    # Create simple diagonal blocks for easy verification
    a = np.array([[1, 0], [0, 1]])
    b = np.array([[2, 0], [0, 2]])
    c = np.array([[3, 0], [0, 3]])
    d = np.array([[4, 0], [0, 4]])

    A = NestedBlockCirculantMatrix([a, b, c, d])
    full_A = np.array(A)

    # Expected structure from docstring:
    # ( a  b | c  d )
    # ( b  a | d  c )
    # ( ----------- )
    # ( c  d | a  b )
    # ( d  c | b  a )
    expected = np.block([
        [a, b, c, d],
        [b, a, d, c],
        [c, d, a, b],
        [d, c, b, a]
    ])

    assert np.allclose(full_A, expected)


def test_nested_block_circulant_lu_decompose():
    """Test LU decomposition works with NestedBlockCirculantMatrix"""
    A = NestedBlockCirculantMatrix([
        2*np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
        np.eye(2) + _rand(2, 2),
    ])
    full_A = np.array(A)

    # LU decompose the converted BlockCirculantMatrix
    bc = A._block_circulant_matrix
    lu_bc = lu_decompose(bc)

    b = RNG.normal(size=(A.shape[0],)) + 1j*RNG.normal(size=(A.shape[0],))
    x_lu = lu_bc.solve(b)
    x_direct = np.linalg.solve(full_A, b)

    assert np.allclose(x_lu, x_direct)


def test_nested_block_circulant_matrix_conversion_is_freed_with_the_matrix():
    # An lru_cache on to_BlockCirculantMatrix used to keep every converted matrix alive
    A = NestedBlockCirculantMatrix([_rand(2, 2) for _ in range(4)])
    bc = A._block_circulant_matrix
    ref_A, ref_bc = weakref.ref(A), weakref.ref(bc)
    del A, bc
    gc.collect()
    assert ref_A() is None and ref_bc() is None


def test_lu_decomposed_matrices_know_their_array_library():
    from capytaine.tools.array_backend import array_namespace
    from capytaine.tools.block_circulant_matrices import BlockDiagonalMatrix
    blocks = [2*np.eye(2) + _rand(2, 2) for _ in range(4)]
    matrices = [
        np.eye(3) + _rand(3, 3),
        BlockDiagonalMatrix(blocks),
        BlockCirculantMatrix(blocks),
        NestedBlockCirculantMatrix(blocks),
    ]
    for A in matrices:
        lu = lu_decompose(A)
        assert lu.__array_namespace__() is array_namespace(np.zeros(1))
        assert lu.device == "cpu"
        assert array_namespace(lu) is array_namespace(np.zeros(1))


def test_block_diagonal_matrix_to_array():
    from capytaine.tools.block_circulant_matrices import BlockDiagonalMatrix
    blocks = [_rand(2, 2) for _ in range(3)]
    A = BlockDiagonalMatrix(blocks)
    expected = np.zeros((6, 6), dtype=complex)
    for i, b in enumerate(blocks):
        expected[2*i:2*i+2, 2*i:2*i+2] = b
    assert np.allclose(np.array(A), expected)
    b = _rand(6)
    assert np.allclose(A.solve(b), np.linalg.solve(expected, b))


def test_blocks_can_be_given_as_an_array():
    blocks = _rand(4, 2, 2)
    assert np.allclose(np.array(BlockCirculantMatrix(blocks)), np.array(BlockCirculantMatrix(list(blocks))))


def test_blocks_given_as_an_array_are_not_copied():
    blocks = _rand(5, 2, 2)
    A = BlockCirculantMatrix(blocks)
    assert A.nb_blocks == 5
    assert np.shares_memory(A.blocks[2], blocks)
    assert np.shares_memory(A.blocks[1:][0], blocks)
    assert A.blocks.as_array(np) is blocks
    assert len(list(A.blocks)) == 5


def test_block_diagonal_matrix_solve_with_several_right_hand_sides():
    from capytaine.tools.block_circulant_matrices import BlockDiagonalMatrix
    blocks = [2*np.eye(2) + _rand(2, 2) for _ in range(3)]
    A = BlockDiagonalMatrix(blocks)
    b = _rand(6, 4)
    assert np.allclose(A.solve(b), np.linalg.solve(np.array(A), b))


def test_block_diagonal_matrix_to_array_keeps_the_precision():
    from capytaine.tools.block_circulant_matrices import BlockDiagonalMatrix
    A = BlockDiagonalMatrix([_rand(2, 2).astype(np.complex64) for _ in range(3)])
    assert np.array(A).dtype == np.complex64
    assert np.array(A, dtype=np.complex128).dtype == np.complex128


def _well_conditioned_blocks(nb_blocks, n=3):
    return [(2*nb_blocks if i == 0 else 1)*np.eye(n) + _rand(n, n) for i in range(nb_blocks)]


# One case for each code path of `block_diagonalize` and `_fourier_shape`:
# explicit formula for 2 blocks or FFT, 2x2 nested blocks or outer FFT only.
CODE_PATHS = [
    (BlockCirculantMatrix, 2),
    (BlockCirculantMatrix, 3),
    (NestedBlockCirculantMatrix, 4),
    (NestedBlockCirculantMatrix, 6),
]


@pytest.mark.parametrize("matrix_type,nb_blocks", CODE_PATHS)
def test_solvers_of_block_circulant_matrices(matrix_type, nb_blocks):
    from capytaine.bem.linear_solvers import solve_gmres
    A = matrix_type(_well_conditioned_blocks(nb_blocks))
    full_A = np.array(A)
    b = _rand(A.shape[0])
    expected = np.linalg.solve(full_A, b)
    assert np.allclose(A.solve(b), expected)
    assert np.allclose(lu_decompose(A).solve(b), expected)
    assert np.allclose(solve_gmres(A, b), expected, atol=1e-5)


def test_fourier_transform_of_rhs_is_consistent_with_block_diagonalization():
    """Check that `_block_fft` with `_fourier_shape` maps A @ x on block_diagonalize(A) @ x_hat."""
    from capytaine.tools.block_circulant_matrices import _block_fft
    for matrix_type, nb_blocks in CODE_PATHS:
        A = matrix_type(_well_conditioned_blocks(nb_blocks))
        x = _rand(A.shape[0])
        assert np.allclose(
            _block_fft(A @ x, A._fourier_shape),
            np.array(A.block_diagonalize()) @ _block_fft(x, A._fourier_shape),
        )


@pytest.mark.parametrize("nb_blocks", [2, 3])
@pytest.mark.parametrize("real_dtype,complex_dtype", [(np.float64, np.complex128), (np.float32, np.complex64)])
def test_block_diagonalize_real_blocks(nb_blocks, real_dtype, complex_dtype):
    blocks = [((2*nb_blocks if i == 0 else 1)*np.eye(3) + RNG.normal(size=(3, 3))).astype(real_dtype)
              for i in range(nb_blocks)]
    A = BlockCirculantMatrix(blocks)
    diag_A = A.block_diagonalize()
    assert diag_A.dtype in (real_dtype, complex_dtype)
    b = RNG.normal(size=A.shape[0])
    assert np.allclose(A.solve(b), np.linalg.solve(np.array(A), b), rtol=1e-4)


def test_gmres_on_block_diagonal_matrix_dispatches_on_the_type_of_the_blocks():
    from capytaine.bem.linear_solvers import solve_gmres
    from capytaine.tools.block_circulant_matrices import BlockDiagonalMatrix
    A = BlockDiagonalMatrix([10*np.eye(3) + _rand(3, 3) for _ in range(4)])
    for b in [_rand(12), np.concatenate([np.zeros(6), _rand(6)])]:
        assert np.allclose(solve_gmres(A, b), np.linalg.solve(np.array(A), b), atol=1e-5)

    # A subclass of np.ndarray standing for the arrays of another array library.
    # (It stays registered after the test, as singledispatch has no unregister.)
    class OtherArray(np.ndarray):
        pass

    calls = []
    @solve_gmres.register(OtherArray)
    def _(A, b):
        calls.append(A.shape)
        return np.linalg.solve(np.asarray(A), b)

    A_other = BlockDiagonalMatrix(np.stack([np.asarray(bl) for bl in A.blocks]).view(OtherArray))
    b = _rand(12)
    assert np.allclose(solve_gmres(A_other, b), np.linalg.solve(np.array(A), b))
    assert calls == [(3, 3)]*4


@pytest.mark.parametrize("matrix_type,nb_blocks", CODE_PATHS)
def test_gmres_keeps_single_precision(matrix_type, nb_blocks):
    from capytaine.bem.linear_solvers import solve_gmres
    A = matrix_type([bl.astype(np.complex64) for bl in _well_conditioned_blocks(nb_blocks)])
    b = _rand(A.shape[0]).astype(np.complex64)
    x = solve_gmres(A, b)
    assert x.dtype == np.complex64
    assert np.allclose(x, np.linalg.solve(np.array(A), b), rtol=1e-3, atol=1e-4)


@pytest.mark.parametrize("matrix_type,nb_blocks", CODE_PATHS)
@pytest.mark.parametrize("real_blocks", [True, False])
def test_solve_with_array_api_strict(matrix_type, nb_blocks, real_blocks):
    xps = pytest.importorskip("array_api_strict")  # Not available on Python 3.8
    blocks = _well_conditioned_blocks(nb_blocks)
    if real_blocks:
        blocks = [bl.real for bl in blocks]
    A = matrix_type(blocks)
    A_strict = matrix_type([xps.asarray(bl) for bl in blocks])
    b = RNG.normal(size=A.shape[0])  # Real right-hand side
    expected = np.linalg.solve(np.array(A), b)
    x = A_strict.solve(xps.asarray(b))
    assert x.__array_namespace__() is xps
    assert np.allclose(np.asarray(x), expected)
    assert np.allclose(np.asarray(A_strict @ xps.asarray(b)), np.array(A) @ b)
