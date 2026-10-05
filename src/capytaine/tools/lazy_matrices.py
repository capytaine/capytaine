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
"""Lazy matrix where the rows are computed and stored on demand when the matrix-vector product is requested."""

from typing import Callable
import numpy as np

from capytaine.tools.array_backend import array_namespace, is_array, to_numpy


def slices(start, stop, chunk_size):
    """Generator returning slices covering the range from start to stop with `chunk_size` elements per slice.

    >>> list(slices(0, 50, 10))
    [slice(0, 15, None),
     slice(15, 30, None),
     slice(30, 45, None),
     slice(45, 50, None)]
    """
    i = start
    while i < stop:
        batch = slice(i, min(i+chunk_size, stop))
        yield batch
        i = i + chunk_size


class LazyMatrix:
    def __init__(self, row_constructor, shape, *, chunk_size=10, dtype=float, namespace=None, device="cpu"):
        """
        A matrix (2D array) that is never fully stored in memory, but instead recomputed from a `row_constructor` method when required.

        Parameters
        ----------
        row_constructor: callable
            Function returning an array containing a few rows of the matrix.
            We assume that row_constructor(slice(n, n+m)) returns an array of shape (m, d) corresponding to the m rows of indices between n and n+m.
            The dtype, array library and device of the output of row_constructor should match the ones given below.
        shape: 2-ple of int
            The shape of the matrix.
        chunk_size: int
            The number of row requested to row_constructor at each call.
        dtype: dtype
            The type of data contained in the matrix.
        namespace: array namespace, optional
            The array library (supporting the array API standard) of the rows (default: NumPy).
        device: device, optional
            The device of the rows (default: "cpu").
        """
        self.row_constructor: Callable[range, np.ndarray] = row_constructor
        self.shape = shape
        self.chunk_size = chunk_size
        self.dtype = dtype
        self.ndim = 2  # Other shapes not implemented
        self._slices = list(slices(0, self.shape[0], self.chunk_size))
        self._namespace = array_namespace(np.empty(0)) if namespace is None else namespace
        self.device = device

    def __array_namespace__(self, *, api_version=None):
        return self._namespace

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        if dtype is None:
            dtype = self.dtype
        rows = [to_numpy(self.row_constructor(sl)) for sl in self._slices]
        return np.concatenate(rows).astype(dtype)

    def __matmul__(self, other):
        if is_array(other) and other.ndim == 1 and other.shape[0] == self.shape[1]:
            # Only matrix-vector product is actually implemented
            # Compute `chunk_size` rows and multiply them by `other`
            output_chunks = [self.row_constructor(sl) @ other for sl in self._slices]
            return self._namespace.concat(output_chunks)
        else:
            return NotImplemented
            # Usually fallback on building the full matrix with __array__ above.
