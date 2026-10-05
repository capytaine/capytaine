# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Helpers to write code that works with any array library implementing the array API standard."""

import numpy as np
from array_api_compat import array_namespace, device, is_numpy_array, to_device

from capytaine.tools.symbolic_multiplication import SymbolicMultiplication

__all__ = ["array_namespace", "is_array", "complex_dtype", "asarray_like", "to_numpy"]


def is_array(x):
    """True if x is an array of any library supported by array-api-compat."""
    try:
        array_namespace(x)
        return True
    except TypeError:
        return False


def complex_dtype(xp, dtype):
    """Complex dtype with the same precision as `dtype`."""
    return xp.complex64 if dtype in (xp.float32, xp.complex64) else xp.complex128


def asarray_like(x, like, *, complex_=False):
    """Convert `x` to the array library, device and (complex) precision of `like`."""
    if isinstance(x, SymbolicMultiplication):
        return SymbolicMultiplication(x.symbol, asarray_like(x.value, like, complex_=complex_))
    xp = array_namespace(like)
    dtype = complex_dtype(xp, like.dtype) if complex_ else like.dtype
    return xp.asarray(x, dtype=dtype, device=device(like))


def to_numpy(x):
    """Convert any array (possibly on GPU) to a NumPy array."""
    if x is None or is_numpy_array(x):
        return x
    if isinstance(x, SymbolicMultiplication):
        return SymbolicMultiplication(x.symbol, to_numpy(x.value))
    try:
        return np.from_dlpack(x)                  # Arrays in host memory, any library
    except (BufferError, RuntimeError, TypeError):
        return np.asarray(to_device(x, "cpu"))    # e.g. torch tensors on GPU
