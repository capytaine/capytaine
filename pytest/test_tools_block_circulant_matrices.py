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

    # Test that to_BlockCirculantMatrix produces equivalent matrix
    bc = A.to_BlockCirculantMatrix()
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

    # Test that to_BlockCirculantMatrix produces equivalent matrix
    bc = A.to_BlockCirculantMatrix()
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

    # Test that to_BlockCirculantMatrix produces equivalent matrix
    bc = A.to_BlockCirculantMatrix()
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
    bc = A.to_BlockCirculantMatrix()
    lu_bc = lu_decompose(bc)

    b = RNG.normal(size=(A.shape[0],)) + 1j*RNG.normal(size=(A.shape[0],))
    x_lu = lu_bc.solve(b)
    x_direct = np.linalg.solve(full_A, b)

    assert np.allclose(x_lu, x_direct)


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
