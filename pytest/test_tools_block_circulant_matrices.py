import numpy as np
from capytaine.tools.block_circulant_matrices import (
        BlockCirculantMatrix, lu_decompose,
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
