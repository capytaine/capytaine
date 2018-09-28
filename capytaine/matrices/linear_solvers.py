#!/usr/bin/env python
# coding: utf-8

import logging
from functools import lru_cache

import numpy as np
from scipy import linalg as sl
from scipy.sparse import linalg as ssl

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import BlockSymmetricToeplitzMatrix, BlockSymmetricCirculantMatrix

LOG = logging.getLogger(__name__)


def solve_directly(A, b):
    if isinstance(A, BlockSymmetricCirculantMatrix):
        LOG.debug("\tSolve linear system %s", A)
        AA = np.array([block.full_matrix() if not isinstance(block, np.ndarray) else block
                       for block in A.first_block_line])
        AAt = np.fft.fft(AA, axis=0)
        bt = np.fft.fft(np.reshape(b, AAt.shape[:2]), axis=0)
        xt = solve_directly(AAt, bt)
        x = np.fft.ifft(xt, axis=0).reshape(b.shape)
        return x

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
    x, info = ssl.gmres(A.full_matrix(), b)
    if info != 0:
        LOG.warning("No convergence of the GMRES")
    return x

