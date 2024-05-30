"""The linear solvers used in Capytaine.

They are based on numpy solvers with a thin layer for the handling of Hierarchical Toeplitz matrices.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

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
            fft_of_result = np.linalg.solve(blocks_of_diagonalization, fft_of_rhs[..., np.newaxis])[..., 0]
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

    if LOG.isEnabledFor(logging.INFO):
        counter = Counter()
        x, info = ssl.gmres(A, b, atol=1e-6, callback=counter)
        LOG.info(f"End of GMRES after {counter.nb_iter} iterations.")

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


# PRECONDITIONED SOLVER

def _block_Jacobi_coarse_corr(A, b, x0, R, RA, AcLU, DLU, diag_shapes, n):
    """
    Performs a single step of the block-Jacobi method with coarse correction.
    Can be used as a preconditioner for matrix A.

    Parameters
    ----------
    A: BlockMatrix
        System matrix
    b: array
        System right hand side vector
    x0: array
        Initial guess of the solution
    R: array
        Coarse space restriction matrix
    RA: array
        Precomputed product of R and A
    AcLU: list
        LU decomposition data of the coarse system matrix Ac (output of
        scipy.linalg.lu_factor)
    DLU: list
        List of LU decomposition data for the diagonal blocks of A
    diag_shapes: list
        List of shapes of diagonal blocks
    n: integer
        Size of the coarse problem (e.g. number of bodies simulated)

    Returns
    -------
    array
        Action of a step of the method on vector x0
    """
    # x_ps = x after the pre-smoothing step (block-Jacobi)
    x_ps = np.zeros(A.shape[0], dtype=complex)

    # the diagonal blocks of A have already been put to zero in build_matrices
    # they are not needed anymore
    q = b - A@x0
    # loop over diagonal blocks
    for kk in range(n):
        local_slice = slice(sum(diag_shapes[:kk]), sum(diag_shapes[:kk+1]))
        local_rhs = q[local_slice]
        local_sol = sl.lu_solve(DLU[kk], local_rhs, check_finite=False)

        x_ps[local_slice] = local_sol

    r_c = R@b - RA@x_ps #restricted residual
    e_c = sl.lu_solve(AcLU, r_c, check_finite=False)
    # update
    return x_ps + R.T@e_c

def solve_precond_gmres(A_and_precond_data, b):
    """
    Implementation of the preconditioner presented in
    `<https://doi.org/10.1007/978-3-031-50769-4_14>`.
    """
    A, R, RA, AcLU, DLU, diag_shapes, n, PinvA = A_and_precond_data
    N = A.shape[0]

    Pinvb = _block_Jacobi_coarse_corr(A, b, np.zeros(N, dtype=complex), R, RA, AcLU, DLU, diag_shapes, n)

    LOG.debug(f"Solve with GMRES for {A}.")

    if LOG.isEnabledFor(logging.INFO):
        counter = Counter()
        x, info = ssl.gmres(PinvA, Pinvb, atol=1e-6, callback=counter)
        LOG.info(f"End of GMRES after {counter.nb_iter} iterations.")

    else:
        x, info = ssl.gmres(PinvA, Pinvb, atol=1e-6)

    if info > 0:
        raise RuntimeError(f"No convergence of the GMRES after {info} iterations.\n"
                            "This can be due to overlapping panels or irregular frequencies.\n"
                            "In the latter case, using a direct solver can help (https://github.com/mancellin/capytaine/issues/30).")

    return x
