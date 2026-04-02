import numpy as np
from capytaine.tools.block_circulant_matrices import (
        BlockCirculantMatrix, NestedBlockCirculantMatrix, lu_decompose,
        leading_dimensions_at_the_end, ending_dimensions_at_the_beginning
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


def test_permute_dims():
    a = RNG.normal(size=(1, 2, 3, 4, 5))
    assert leading_dimensions_at_the_end(a).shape == (3, 4, 5, 1, 2)
    assert ending_dimensions_at_the_beginning(a).shape == (4, 5, 1, 2, 3)
    assert np.allclose(ending_dimensions_at_the_beginning(
        leading_dimensions_at_the_end(a)
    ), a)
    assert np.allclose(leading_dimensions_at_the_end(
        ending_dimensions_at_the_beginning(a)
    ), a)


def test_deeper_2x2_block_circulant_matrices():
    A = BlockCirculantMatrix([
        _rand(2, 2, 3),
        _rand(2, 2, 3),
    ])
    full_A = np.array(A)
    assert full_A.shape == (4, 4, 3)


def test_2x2_nested_block_circulant_matrices():
    A = BlockCirculantMatrix([
            BlockCirculantMatrix([
                2*np.eye(2) + _rand(2, 2),
                np.zeros((2, 2)) + _rand(2, 2),
            ]),
            BlockCirculantMatrix([
                np.eye(2) + _rand(2, 2),
                np.zeros((2, 2)) + _rand(2, 2),
            ]),
        ])
    full_A = np.array(A)
    assert full_A.shape == (8, 8)
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

def test_nested_2x2_3x3_block_circulant_matrices():
    A = BlockCirculantMatrix([
        BlockCirculantMatrix([
            1*np.eye(2) + _rand(2, 2),
            2*np.eye(2) + _rand(2, 2),
            3*np.eye(2) + _rand(2, 2)
            ]),
        BlockCirculantMatrix([
            np.zeros(2) + _rand(2, 2),
            np.zeros(2) + _rand(2, 2),
            np.zeros(2) + _rand(2, 2)
            ]),
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

def test_nested_3x3_2x2_block_circulant_matrices():
    A = BlockCirculantMatrix([
        BlockCirculantMatrix([
            1*np.eye(2) + _rand(2, 2),
            2*np.eye(2) + _rand(2, 2)
            ]),
        BlockCirculantMatrix([
            np.zeros(2) + _rand(2, 2),
            np.zeros(2) + _rand(2, 2)
            ]),
        BlockCirculantMatrix([
            np.zeros(2) + _rand(2, 2),
            np.zeros(2) + _rand(2, 2)
            ]),
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
