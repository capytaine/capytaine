#!/usr/bin/env python
# coding: utf-8
"""This module implements several classes describing matrices defined by blocks.
These matrices can be nested to recursively define Hierarchical matrices.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.block_toeplitz import (
    BlockToeplitzMatrix, BlockSymmetricToeplitzMatrix,
    BlockCirculantMatrix, EvenBlockSymmetricCirculantMatrix, OddBlockSymmetricCirculantMatrix,
)
from capytaine.matrices.builders import (
    cut_matrix, random_block_matrix,
    full_like, zeros_like, ones_like, identity_like,
)
from capytaine.matrices.low_rank import LowRankMatrix
