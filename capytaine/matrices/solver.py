#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import BlockSymmetricToeplitzMatrix, BlockSymmetricCirculantMatrix

LOG = logging.getLogger(__name__)


def solve(A, b):
    if isinstance(A, BlockSymmetricCirculantMatrix):
        LOG.debug("\tSolve linear system %s BlockCirculantMatrix (block size: %s)",
                  A.nb_blocks, A.block_shape)
        AA = np.array([block.full_matrix() if not isinstance(block, np.ndarray) else block
                       for block in A.first_block_line])
        AAt = np.fft.fft(AA, axis=0)
        bt = np.fft.fft(np.reshape(b, AAt.shape[:2]), axis=0)
        xt = solve(AAt, bt)
        x = np.fft.ifft(xt, axis=0).reshape(b.shape)
        return x

    elif isinstance(A, BlockSymmetricToeplitzMatrix):
        if A.nb_blocks == (2, 2):
            LOG.debug("\tSolve system of 2×2 BlockToeplitzMatrix (block size: %s)", A.block_shape)
            A1, A2 = A.first_block_line
            b1, b2 = b[:len(b)//2], b[len(b)//2:]
            x_plus = solve(A1 + A2, b1 + b2)
            x_minus = solve(A1 - A2, b1 - b2)
            return np.concatenate([x_plus + x_minus, x_plus - x_minus])/2
        else:
            # Not implemented
            return solve(A.full_matrix(), b)

    elif isinstance(A, BlockMatrix):
        return solve(A.full_matrix(), b)

    elif isinstance(A, np.ndarray):
        LOG.debug(f"\tSolve linear system (size: {A.shape}) with numpy.")
        return np.linalg.solve(A, b)

    else:
        raise ValueError(f"Unrecognized type of matrix to solve: {A}")

