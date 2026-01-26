"""This module implements several classes describing matrices defined by blocks.
These matrices can be nested to recursively define Hierarchical matrices.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/fakeblocks>

from fakeblocks.matrices.block import BlockMatrix
from fakeblocks.matrices.block_toeplitz import (
    BlockToeplitzMatrix, BlockSymmetricToeplitzMatrix,
    BlockCirculantMatrix, EvenBlockSymmetricCirculantMatrix, OddBlockSymmetricCirculantMatrix,
)
from fakeblocks.matrices.builders import (
    cut_matrix, random_block_matrix,
    full_like, zeros_like, ones_like, identity_like,
)
from fakeblocks.matrices.low_rank import LowRankMatrix
