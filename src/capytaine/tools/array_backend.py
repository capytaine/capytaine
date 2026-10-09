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

from typing import Any, Protocol, Tuple

import numpy as np
from array_api_compat import array_namespace, device, is_numpy_array, to_device

from capytaine.tools.symbolic_multiplication import SymbolicMultiplication

__all__ = [
    "MatrixLike", "array_namespace", "device", "is_array", "is_numpy_namespace", "complex_dtype", "to_backend_of", "to_numpy",
    "IterableTensor", "split", "leading_dimensions_at_the_end", "ending_dimensions_at_the_beginning",
]


class MatrixLike(Protocol):
    """Any 2D (or batched) matrix: an array of any library supported by array-api-compat,
    or a custom data-sparse matrix such as BlockCirculantMatrix.

    Only used for type annotations. `__array_namespace__` is not required, since
    some arrays (e.g. torch tensors) get it from array-api-compat instead of defining it."""
    @property
    def shape(self) -> Tuple[int, ...]: ...
    @property
    def dtype(self) -> Any: ...
    def __matmul__(self, other: Any) -> Any: ...


def is_array(x):
    """True if x is an array of any library supported by array-api-compat."""
    try:
        array_namespace(x)
        return True
    except TypeError:
        return False


def is_numpy_namespace(xp):
    """True if `xp` is the NumPy namespace.
    (`array_api_compat.is_numpy_namespace` is only available since array-api-compat 1.9.)"""
    return xp is array_namespace(np.empty(0))


def complex_dtype(xp, dtype):
    """Complex dtype with the same precision as `dtype`."""
    return xp.complex64 if dtype in (xp.float32, xp.complex64) else xp.complex128


def to_backend_of(x, like, *, complex_=False):
    """Convert `x` to the array library, device and (complex) precision of `like`."""
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


def split(x, n):
    """Split an array in `n` parts of equal sizes along its first axis (replaces `np.split`, that is not in the array API standard)."""
    size = x.shape[0] // n
    return [x[i*size:(i+1)*size, ...] for i in range(n)]


def leading_dimensions_at_the_end(a):
    """Transform an array of shape (n, m, ...) into (..., n, m).
    Invert of `ending_dimensions_at_the_beginning`"""
    xp = array_namespace(a)
    return xp.permute_dims(a, (*range(2, a.ndim), 0, 1))


def ending_dimensions_at_the_beginning(a):
    """Transform an array of shape (..., n, m) into (n, m, ...).
    Invert of `leading_dimensions_at_the_end`"""
    xp = array_namespace(a)
    return xp.permute_dims(a, (a.ndim - 2, a.ndim - 1, *range(a.ndim - 2)))


class IterableTensor:
    """Sequence of the sub-arrays along the first axis of an array, or of a list of arrays.

    It can be built from a list of arrays of the same shape or from a single array
    of shape (n, ...), which is then kept as it is (no copy, e.g. to apply a FFT across the first axis)
    and the items are the views `array[i, ...]`.

    Arrays of the array API standard can not always be iterated, sliced, unpacked or measured with `len`
    (unlike the ones of NumPy, PyTorch or JAX), and an index has to cover all the dimensions of the array.
    This class gives the usual behavior of a sequence in all cases.
    """
    def __init__(self, tensor):
        self._array = tensor if is_array(tensor) else None
        self._list = None if is_array(tensor) else list(tensor)

    def __len__(self):
        return self._array.shape[0] if self._array is not None else len(self._list)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        return self._array[index, ...] if self._array is not None else self._list[index]

    def __iter__(self):
        return (self[i] for i in range(len(self)))

    def as_array(self, xp):
        """All the items in a single array of shape (n, ...)"""
        return self._array if self._array is not None else xp.stack(self._list)
