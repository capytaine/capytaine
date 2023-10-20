"""Prony decomposition: tool to approximate a function as a sum of exponentials.
Used in particular in the finite depth Green function.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
from numpy.polynomial import polynomial
from scipy.optimize import curve_fit
from scipy.linalg import toeplitz

LOG = logging.getLogger(__name__)


def exponential_decomposition(X, F, m):
    """Use Prony's method to approximate the sampled real function F=f(X) as a sum of m
    exponential functions x → Σ a_i exp(lamda_i x).

    Parameters
    ----------
    X: 1D array
        sampling points.
    F: 1D array (same size as X)
        values of the function to approximate at the points of x.
    m: integer
        number of exponential functions

    Return
    ------
    a: 1D array (size m)
        coefficients of the exponentials
    lamda: 1D array (size m)
        growth rate of the exponentials
    """
    assert X.shape == F.shape

    # Compute the coefficients of the polynomials of Prony's method
    A = toeplitz(c=F[m-1:-1], r=F[:m][::-1])
    P, *_ = np.linalg.lstsq(A, F[m:], rcond=None)

    # Build and solve polynomial function
    coeffs = np.ones(m+1)
    # coeffs[:m] = -P[::-1]
    for i in range(m):
        coeffs[m-i-1] = -P[i]
    roots = polynomial.polyroots(coeffs)

    # Discard values where log is undefined
    roots = roots[np.logical_or(np.imag(roots) != 0.0, np.real(roots) >= 0.0)]

    # Deduce lamda and keep only interesting values
    lamda = np.real(np.log(roots)/(X[1] - X[0]))
    lamda = np.unique(lamda)
    lamda = lamda[np.logical_and(-20.0 < lamda, lamda < 0.0)]

    # Fit the values of 'a' on the curve
    def f(x, *ar):
        ar = np.asarray(ar)[:, np.newaxis]
        la = lamda[:, np.newaxis]
        return np.sum(ar * np.exp(la * x), axis=0)
    a, *_ = curve_fit(f, X, F, p0=np.zeros(lamda.shape))

    return a, lamda


def error_exponential_decomposition(X, F, a, lamda):
    """Compare exponential decomposition defined by the coefficients a and lamda to the reference
    values in F.

    Parameters
    ----------
    X: 1D array
        sampling points
    F: 1D array (same size as X)
        reference values
    a: 1D array
        coefficients of the exponentials
    lamda: 1D array (same size as a)
        growth rate of the exponentials

    Returns
    -------
    error: float
        mean square error of the decomposition
    """
    a = np.asarray(a)[:, np.newaxis]
    lamda = np.asarray(lamda)[:, np.newaxis]

    def f(x):
        return np.sum(a * np.exp(lamda*x), axis=0)

    return np.square(f(X) - F).mean()
