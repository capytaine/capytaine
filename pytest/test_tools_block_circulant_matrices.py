import numpy as np
from capytaine.tools.block_circulant_matrices import BlockCirculantMatrix, lu_decompose

RNG = np.random.default_rng(seed=0)

def test_2x2_block_circulant_matrices():
    A = BlockCirculantMatrix([
        2*np.eye(2) + RNG.normal(size=(2, 2)),
        np.eye(2) + RNG.normal(size=(2, 2)),
    ])
    full_A = np.array(A)
    assert full_A.shape == (4, 4)

    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )

def test_2x2_nested_block_circulant_matrices():
    A = BlockCirculantMatrix([
            BlockCirculantMatrix([
                2*np.eye(2) + RNG.normal(size=(2, 2)),
                np.zeros((2, 2)) + RNG.normal(size=(2, 2)),
            ]),
            BlockCirculantMatrix([
                np.eye(2) + RNG.normal(size=(2, 2)),
                np.zeros((2, 2)) + RNG.normal(size=(2, 2)),
            ]),
        ])
    full_A = np.array(A)
    assert full_A.shape == (8, 8)
    b = RNG.normal(size=(A.shape[0],))
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
        np.eye(2) + RNG.normal(size=(2, 2)),
        2*np.eye(2) + RNG.normal(size=(2, 2)),
        3*np.eye(2) + RNG.normal(size=(2, 2))
    ])
    full_A = np.array(A)
    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
        lu_decompose(A).solve(b),
        np.linalg.solve(full_A, b)
    )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )

def test_nested_3x3_block_circulant_matrices():
    A = BlockCirculantMatrix([
        BlockCirculantMatrix([
            1*np.eye(2) + RNG.normal(size=(2, 2)),
            2*np.eye(2) + RNG.normal(size=(2, 2)),
            3*np.eye(2) + RNG.normal(size=(2, 2))
            ]),
        BlockCirculantMatrix([
            np.zeros(2) + RNG.normal(size=(2, 2)),
            np.zeros(2) + RNG.normal(size=(2, 2)),
            np.zeros(2) + RNG.normal(size=(2, 2))
            ]),
        ])
    full_A = np.array(A)
    b = RNG.normal(size=(A.shape[0],))
    assert np.allclose(
            lu_decompose(A).solve(b),
            np.linalg.solve(full_A, b)
            )
    assert np.allclose(
        lu_decompose(A).solve(b),
        A.solve(b)
    )
