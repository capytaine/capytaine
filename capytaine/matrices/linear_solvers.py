#!/usr/bin/env python
# coding: utf-8

import logging
from functools import lru_cache

import numpy as np
from scipy import linalg as sl
from scipy.sparse import linalg as ssl

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import (
    BlockSymmetricToeplitzMatrix,
    _AbstractBlockSymmetricCirculantMatrix
)

LOG = logging.getLogger(__name__)


def solve_directly(A, b):
    assert isinstance(b, np.ndarray) and A.ndim == b.ndim+1 and A.shape[-2] == b.shape[-1]
    if isinstance(A, _AbstractBlockSymmetricCirculantMatrix):
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
            A1, A2 = A.first_block_line
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


@lru_cache(maxsize=1)
def lu_decomp(A):
    LOG.debug(f"Compute LU decomposition of {A}.")
    return sl.lu_factor(A.full_matrix())


def solve_storing_lu(A, b):
    LOG.debug(f"Solve with LU decomposition of {A}.")
    return sl.lu_solve(lu_decomp(A), b)


def solve_gmres(A, b):
    LOG.debug(f"Solve with GMRES for {A}.")
    x, info = ssl.gmres(A, b, atol=1e-6)
    if info != 0:
        LOG.warning("No convergence of the GMRES")
    return x

