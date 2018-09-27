#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np

from capytaine.matrices.block_matrices import BlockMatrix
from capytaine.matrices.block_toeplitz_matrices import BlockSymmetricToeplitzMatrix, BlockSymmetricCirculantMatrix

LOG = logging.getLogger(__name__)


def solve(A, b):

    if isinstance(A, BlockSymmetricToeplitzMatrix):
        return solve(A.full_matrix(), b)

    elif isinstance(A, BlockSymmetricToeplitzMatrix):
        return solve(A.full_matrix(), b)

    elif isinstance(A, BlockMatrix):
        return solve(A.full_matrix(), b)

    elif isinstance(A, np.ndarray):
        LOG.debug(f"\tSolve linear system (size: {A.shape}) with numpy.")
        return np.linalg.solve(A, b)

    else:
        raise ValueError(f"Unrecognized type of {A} in solve")

