import logging
from abc import ABC, abstractmethod
from functools import singledispatch
from typing import Any

import numpy as np
from numpy.typing import NDArray
import scipy.linalg as sl
import scipy.sparse.linalg as ssl

from capytaine.tools.array_backend import  array_namespace, device


LOG = logging.getLogger(__name__)


############################### LU DECOMPOSITION ###############################

# ABSTRACT INTERFACE FOR LU DECOMPOSED MATRICES

class AbstractLUDecomposedMatrix(ABC):
    """Base class of the LU decompositions of matrices.

    New matrix types (e.g. matrices stored with another array library) can be
    supported by the linear solvers of Capytaine by registering a function
    returning a subclass of this class with :func:`lu_decompose`.

    As for an array of the array API standard, the subclasses define the method
    `__array_namespace__` and the attribute (or property) `device`, returning
    the array library and the device of the matrix that has been decomposed.
    """
    shape: tuple
    dtype: object
    device: Any

    @abstractmethod
    def __array_namespace__(self, *, api_version=None):
        """Array library of the matrix that has been decomposed, as for an array of the array API standard."""

    @abstractmethod
    def solve(self, b):
        """Solve the linear system with the decomposed matrix as left-hand side and `b` as right-hand side."""


def has_been_lu_decomposed(A):
    return isinstance(A, AbstractLUDecomposedMatrix)


@singledispatch
def lu_decompose(A, *, overwrite_a: bool = False) -> AbstractLUDecomposedMatrix:
    """Compute the LU decomposition of `A`.

    This is a :func:`functools.singledispatch` function: the implementation for
    a new type of matrix can be added with ``lu_decompose.register(MyMatrix, my_function)``,
    where ``my_function(A, *, overwrite_a=False)`` returns an :class:`AbstractLUDecomposedMatrix`.
    """
    raise NotImplementedError(f"No LU decomposition registered for {type(A)}")


# IMPLEMENTATION FOR NUMPY ARRAYS

class LUDecomposedMatrix(AbstractLUDecomposedMatrix):
    def __init__(self, A: NDArray, *, overwrite_a : bool = False):
        LOG.debug("LU decomp of %s of shape %s",
                  A.__class__.__name__, A.shape)
        self._lu_decomp = sl.lu_factor(A, overwrite_a=overwrite_a)
        self.shape = A.shape
        self.dtype = A.dtype
        self._array_backend = array_namespace(A)
        self.device = device(A)

    def __array_namespace__(self, *, api_version=None):
        return self._array_backend

    def solve(self, b: np.ndarray) -> np.ndarray:
        LOG.debug("Called solve on %s of shape %s",
                  self.__class__.__name__, self.shape)
        return sl.lu_solve(self._lu_decomp, b)


@lu_decompose.register(np.ndarray)
def _lu_decompose_ndarray(A, *, overwrite_a: bool = False):
    return LUDecomposedMatrix(A, overwrite_a=overwrite_a)

#################################### GMRES #####################################

@singledispatch
def solve_gmres(A, b):
    """Solve `A x = b` with GMRES, where `b` is a vector.

    This is a :func:`functools.singledispatch` function: the implementation is chosen with the type of `A`.
    GMRES for a new type of matrices (e.g. arrays of another array library) can be supported by
    registering an implementation with ``@solve_gmres.register(MatrixType)``.
    """
    raise NotImplementedError(
        f"The 'gmres' linear solver has not been implemented for matrix of {type(A)}"
    )


class _Counter:
    def __init__(self):
        self.nb_iter = 0

    def __call__(self, *args, **kwargs):
        self.nb_iter += 1


@solve_gmres.register(np.ndarray)
def _solve_gmres_scipy(A, b):
    LOG.debug(f"Solve with GMRES for {A}.")

    if LOG.isEnabledFor(logging.INFO):
        counter = _Counter()
        x, info = ssl.gmres(A, b, atol=1e-6, callback=counter, callback_type="legacy")
        LOG.info(f"End of GMRES after {counter.nb_iter} iterations.")

    else:
        x, info = ssl.gmres(A, b, atol=1e-6)

    if info > 0:
        raise RuntimeError(f"No convergence of the GMRES after {info} iterations.\n"
                            "This can be due to overlapping panels or irregular frequencies.")

    return x
