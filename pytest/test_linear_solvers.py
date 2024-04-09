import pytest
import numpy as np
from numpy.random import default_rng
import capytaine as cpt

from capytaine.matrices.block_toeplitz import BlockSymmetricToeplitzMatrix, BlockCirculantMatrix, BlockToeplitzMatrix
from capytaine.matrices.builders import random_block_matrix
from capytaine.matrices.linear_solvers import solve_directly, LUSolverWithCache, solve_gmres

RNG = default_rng(seed=0)

#######################################################################
#                            Full problems                            #
#######################################################################

@pytest.fixture
def solved_full_problem():
    A = RNG.random((5, 5))
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_full_problem(solved_full_problem):
    A, x_ref, b = solved_full_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_solve_with_lu_full_problem(solved_full_problem):
    A, x_ref, b = solved_full_problem
    linear_solver = LUSolverWithCache()
    x = linear_solver.solve(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)
    x = linear_solver.solve(A, b)  # Should use cache
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_gmres_full_problem(solved_full_problem):
    A, x_ref, b = solved_full_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

#######################################################################
#           2x2 block symmetric matrices (reflection mesh)            #
#######################################################################

@pytest.fixture
def solved_2x2_block_symmetric_toeplitz_problem():
    A = BlockSymmetricToeplitzMatrix([
        [RNG.random((4, 4)) for _ in range(2)]
    ])
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_2x2_block_symmetric_toeplitz_problem(solved_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_2x2_block_symmetric_toeplitz_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_solve_with_lu_2x2_block_symmetric_toeplitz_problem(solved_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_2x2_block_symmetric_toeplitz_problem
    linear_solver = LUSolverWithCache()
    x = linear_solver.solve(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)
    x = linear_solver.solve(A, b)  # Should use cache
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_gmres_2x2_block_symmetric_toeplitz_problem(solved_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_2x2_block_symmetric_toeplitz_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

#######################################################################
#       Nested 2x2 block symmetric matrices (two reflections)         #
#######################################################################

@pytest.fixture
def solved_nested_2x2_block_symmetric_toeplitz_problem():
    A = BlockSymmetricToeplitzMatrix([
        [BlockSymmetricToeplitzMatrix([
            [RNG.random((4, 4)) for _ in range(2)]
            ]) for _ in range(2)]
        ])
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_nested_2x2_block_symmetric_toeplitz_problem(solved_nested_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_nested_2x2_block_symmetric_toeplitz_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_solve_with_lu_nested_2x2_block_symmetric_toeplitz_problem(solved_nested_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_nested_2x2_block_symmetric_toeplitz_problem
    linear_solver = LUSolverWithCache()
    x = linear_solver.solve(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)
    x = linear_solver.solve(A, b)  # Should use cache
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_gmres_2x2_nested_block_symmetric_toeplitz_problem(solved_nested_2x2_block_symmetric_toeplitz_problem):
    A, x_ref, b = solved_nested_2x2_block_symmetric_toeplitz_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

#######################################################################
#                           Block circulant                           #
#######################################################################

@pytest.fixture
def solved_block_circulant_problem():
    A = BlockCirculantMatrix([
        [RNG.random((3, 3)) for _ in range(6)]
    ])
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_block_circulant_problem(solved_block_circulant_problem):
    A, x_ref, b = solved_block_circulant_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_solve_with_lu_block_circulant_problem(solved_block_circulant_problem):
    A, x_ref, b = solved_block_circulant_problem
    linear_solver = LUSolverWithCache()
    with pytest.raises(NotImplementedError):
        linear_solver.solve(A, b)

def test_gmres_block_circulant_problem(solved_block_circulant_problem):
    A, x_ref, b = solved_block_circulant_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)


#######################################################################
#                           Block Toeplitz                            #
#######################################################################
@pytest.fixture
def solved_block_toeplitz_problem():
    A = BlockToeplitzMatrix([
        [RNG.random((3, 3)) for _ in range(7)]
    ])
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_block_toeplitz_problem(solved_block_toeplitz_problem):
    A, x_ref, b = solved_block_toeplitz_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_gmres_block_toeplitz_problem(solved_block_toeplitz_problem):
    A, x_ref, b = solved_block_toeplitz_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)


#######################################################################
#                            Nested blocks                            #
#######################################################################
@pytest.fixture
def solved_nested_block_problem():
    A = BlockCirculantMatrix([
        [random_block_matrix([1, 1], [1, 1], rng=RNG) for _ in range(6)]
    ])
    x = RNG.random(A.shape[0])
    return (A, x, A @ x)

def test_solve_directly_nested_block_problem(solved_nested_block_problem):
    A, x_ref, b = solved_nested_block_problem
    x = solve_directly(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)

def test_gmres_nested_block_problem(solved_nested_block_problem):
    A, x_ref, b = solved_nested_block_problem
    x = solve_gmres(A, b)
    assert np.allclose(x, x_ref, rtol=1e-10)
