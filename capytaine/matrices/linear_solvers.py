#!/usr/bin/env python
# coding: utf-8
"""The linear solvers used in Capytaine.

They are based on numpy solvers with a thin layer for the handling of Hierarchical Toeplitz matrices.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from functools import lru_cache

import numpy as np
from scipy import linalg as sl
from scipy.sparse import linalg as ssl

from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.block_toeplitz import BlockSymmetricToeplitzMatrix, BlockCirculantMatrix

LOG = logging.getLogger(__name__)


# DIRECT SOLVER

def solve_directly(A, b):
    assert isinstance(b, np.ndarray) and A.ndim == b.ndim+1 and A.shape[-2] == b.shape[-1]
    if isinstance(A, BlockCirculantMatrix):
        LOG.debug("\tSolve linear system %s", A)
        blocks_of_diagonalization = A.block_diagonalize()
        fft_of_rhs = np.fft.fft(np.reshape(b, (A.nb_blocks[0], A.block_shape[0])), axis=0)
        try:  # Try to run it as vectorized numpy arrays.
            fft_of_result = np.linalg.solve(blocks_of_diagonalization, fft_of_rhs)
        except np.linalg.LinAlgError:  # Or do the same thing with list comprehension.
            fft_of_result = np.array([solve_directly(block, vec) for block, vec in zip(blocks_of_diagonalization, fft_of_rhs)])
        result = np.fft.ifft(fft_of_result, axis=0).reshape((A.shape[1],))
        return result

    elif isinstance(A, BlockSymmetricToeplitzMatrix):
        if A.nb_blocks == (2, 2):
            LOG.debug("\tSolve linear system %s", A)
            A1, A2 = A._stored_blocks[0, :]
            b1, b2 = b[:len(b)//2], b[len(b)//2:]
            x_plus = solve_directly(A1 + A2, b1 + b2)
            x_minus = solve_directly(A1 - A2, b1 - b2)
            return np.concatenate([x_plus + x_minus, x_plus - x_minus])/2
        else:
            # Not implemented
            LOG.debug("\tSolve linear system %s", A)
            return solve_directly(A.full_matrix(), b)

    elif isinstance(A, BlockMatrix):
        LOG.debug("\tSolve linear system %s", A)
        return solve_directly(A.full_matrix(), b)

    elif isinstance(A, np.ndarray):
        LOG.debug(f"\tSolve linear system (size: {A.shape}) with numpy direct solver.")
        return np.linalg.solve(A, b)

    else:
        raise ValueError(f"Unrecognized type of matrix to solve: {A}")


# CACHED LU DECOMPOSITION
class LUSolverWithCache:
    """Solve linear system with the LU decomposition.

    The latest LU decomposition is kept in memory, if a system with the same matrix needs to be solved again, then the decomposition is reused.

    Most of the complexity of this class comes from:
    1. @lru_cache does not work because numpy arrays are not hashable. So a basic cache system has been recoded from scratch.
    2. To be the default solver for the BasicMatrixEngine, the solver needs to support matrices for problems with one or two reflection symmetries.
       Hence, a custom way to cache the LU decomposition of the matrices involved in the direct linear resolution of the symmetric problem.
    """
    def __init__(self):
        self.cached_matrix = None
        self.cached_decomp = None

    def solve(self, A, b):
        return self.solve_with_decomp(self.cached_lu_decomp(A), b)

    def lu_decomp(self, A):
        """Return the LU decomposition of A.
        If A is BlockSymmetricToeplitzMatrix, then return a list of LU decompositions for each block of the block diagonalisation of the matrix.
        """
        if isinstance(A, BlockSymmetricToeplitzMatrix) and A.nb_blocks == (2, 2):
            A1, A2 = A._stored_blocks[0, :]
            return [self.lu_decomp(A1 + A2), self.lu_decomp(A1 - A2)]
        elif isinstance(A, np.ndarray):
            return sl.lu_factor(A)
        else:
            raise NotImplementedError("Cached LU solver is only implemented for dense matrices and 2×2 BlockSymmetricToeplitzMatrix.")

    def cached_lu_decomp(self, A):
        if not(A is self.cached_matrix):
            self.cached_matrix = A
            LOG.debug(f"Computing and caching LU decomposition")
            self.cached_decomp = self.lu_decomp(A)
        else:
            LOG.debug(f"Using cached LU decomposition")
        return self.cached_decomp

    def solve_with_decomp(self, decomp, b):
        """Solve the system using the precomputed LU decomposition.
        TODO: find a better way to differentiate a LU decomposition (returned as tuple by sl.lu_factor)
        and a set of LU decomposition (stored as a list in self.lu_decomp).
        """
        if isinstance(decomp, list):  # The matrix was a BlockSymmetricToeplitzMatrix
            b1, b2 = b[:len(b)//2], b[len(b)//2:]
            x_plus = self.solve_with_decomp(decomp[0], b1 + b2)
            x_minus = self.solve_with_decomp(decomp[1], b1 - b2)
            return np.concatenate([x_plus + x_minus, x_plus - x_minus])/2
        elif isinstance(decomp, tuple):  # The matrix was a np.ndarray
            return sl.lu_solve(decomp, b)
        else:
            raise NotImplementedError("Cached LU solver is only implemented for dense matrices and 2×2 BlockSymmetricToeplitzMatrix.")


# ITERATIVE SOLVER

class Counter:
    def __init__(self):
        self.nb_iter = 0

    def __call__(self, *args, **kwargs):
        self.nb_iter += 1


def solve_gmres(A, b):
    LOG.debug(f"Solve with GMRES for {A}.")

    if LOG.isEnabledFor(logging.DEBUG):
        counter = Counter()
        x, info = ssl.gmres(A, b, atol=1e-6, callback=counter)
        LOG.debug(f"End of GMRES after {counter.nb_iter} iterations.")

    else:
        x, info = ssl.gmres(A, b, atol=1e-6)

    if info > 0:
        raise RuntimeError(f"No convergence of the GMRES after {info} iterations.\n"
                            "This can be due to overlapping panels or irregular frequencies.\n"
                            "In the latter case, using a direct solver can help (https://github.com/mancellin/capytaine/issues/30).")

    return x

def gmres_no_fft(A, b):
    LOG.debug(f"Solve with GMRES for {A} without using FFT.")

    x, info = ssl.gmres(A.no_toeplitz() if isinstance(A, BlockMatrix) else A, b, atol=1e-6)

    if info != 0:
        LOG.warning(f"No convergence of the GMRES. Error code: {info}")

    return x

