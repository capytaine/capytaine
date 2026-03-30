"""Lazy matrix where the rows are computed and stored on demand when the matrix-vector product is requested."""

from typing import Callable
import numpy as np


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
    def __init__(self, row_constructor, shape, *, chunk_size=10, dtype=float):
        """
        A matrix (2D array) that is never fully stored in memory, but instead recomputed from a `row_constructor` method when required.

        Parameters
        ----------
        row_constructor: callable
            Function returning a numpy array containing a few rows of the matrix.
            We assume that row_constructor(slice(n, n+m)) returns a numpy array of shape (m, d) corresponding to the m rows of indices between n and n+m.
            The dtype of the output of row_constructor should match self.dtype.
        shape: 2-ple of int
            The shape of the matrix.
        chunk_size: int
            The number of row requested to row_constructor at each call.
        dtype: numpy.dtype
            The type of data contained in the matrix.
        """
        self.row_constructor: Callable[range, np.ndarray] = row_constructor
        self.shape = shape
        self.chunk_size = chunk_size
        self.dtype = dtype
        self.ndim = 2  # Other shapes not implemented
        self._slices = list(slices(0, self.shape[0], self.chunk_size))

    def __array__(self, dtype=None, copy=True):
        if not copy:
            raise NotImplementedError
        if dtype is None:
            dtype = self.dtype
        rows = [self.row_constructor(sl) for sl in self._slices]
        return np.concatenate(rows).astype(dtype)

    def __matmul__(self, other):
        if isinstance(other, np.ndarray) and other.ndim == 1 and other.shape[0] == self.shape[1]:
            # Only matrix-vector product is actually implemented
            # Compute `chunk_size` rows and multiply them by `other`
            output_chunks = [self.row_constructor(sl) @ other for sl in self._slices]
            return np.concatenate(output_chunks)
        else:
            return NotImplemented
            # Usually fallback on building the full matrix with __array__ above.
