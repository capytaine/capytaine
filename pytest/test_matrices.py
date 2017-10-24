#!/usr/bin/env python
# coding: utf-8

import pytest

from capytaine.Toeplitz_matrices import *


def test_BlockToeplitz():
    A = BlockToeplitzMatrix([np.array([[i]]) for i in range(5)])
    assert A.nb_blocks == 5
    assert A.block_size == 1
    assert np.all(A.full_matrix() == np.array([[0, 1, 2, 3, 4],
                                               [1, 0, 1, 2, 3],
                                               [2, 1, 0, 1, 2],
                                               [3, 2, 1, 0, 1],
                                               [4, 3, 2, 1, 0]])
                  )

    A = BlockToeplitzMatrix([np.zeros((3, 3)) for _ in range(10)])
    assert A.nb_blocks == 10
    assert A.block_size == 3
    assert np.all(A.full_matrix() == np.zeros((30, 30)))
