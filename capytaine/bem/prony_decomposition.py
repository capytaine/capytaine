#!/usr/bin/env python
# coding: utf-8
"""Prony decomposition: tool to approximate a function as a sum of exponentials.
Used in particular in the finite depth Green function.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from functools import lru_cache

import numpy as np
from numpy.polynomial import polynomial
from scipy.optimize import curve_fit
from scipy.linalg import toeplitz

import capytaine.bem.NemohCore as NemohCore

LOG = logging.getLogger(__name__)


@lru_cache(maxsize=128)
def find_best_exponential_decomposition(dimensionless_omega, dimensionless_wavenumber, method='fortran'):
    """Compute the decomposition of a part of the finite depth Green function as a sum of exponential functions.

    Two implementations are available: the legacy Fortran implementation from Nemoh and a newer one written in Python.
    For some still unexplained reasons, the two implementations do not always give the exact same result.
    Until the problem is better understood, the Fortran implementation is the default one, to ensure consistency with Nemoh.
    The Fortran version is also significantly faster...

    Results are cached.

    Parameters
    ----------
    dimensionless_omega: float
        dimensionless angular frequency: :math:`kh \\tanh (kh) = \omega^2 h/g`
    dimensionless_wavenumber: float
        dimensionless wavenumber: :math:`kh`
    method: string, optional
        the implementation that should be used to compute the Prony decomposition

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        the amplitude and growth rates of the exponentials
    """

    LOG.debug(f"\tCompute Prony decomposition in finite depth Green function "
              f"for dimless_omega=%.2e and dimless_wavenumber=%.2e",
              dimensionless_omega, dimensionless_wavenumber)

    if method.lower() == 'python':
        # The function that will be approximated.
        @np.vectorize
        def f(x):
            return NemohCore.initialize_green_wave.ff(x, dimensionless_omega, dimensionless_wavenumber)

        # Try different increasing number of exponentials
        for n_exp in range(4, 31, 2):

            # The coefficients are computed on a resolution of 4*n_exp+1 ...
            X = np.linspace(-0.1, 20.0, 4*n_exp+1)
            a, lamda = exponential_decomposition(X, f(X), n_exp)

            # ... and they are evaluated on a finer discretization.
            X = np.linspace(-0.1, 20.0, 8*n_exp+1)
            if error_exponential_decomposition(X, f(X), a, lamda) < 1e-4:
                break

        else:
            LOG.warning("No suitable exponential decomposition has been found"
                        "for dimless_omega=%.2e and dimless_wavenumber=%.2e",
                        dimensionless_omega, dimensionless_wavenumber)

    elif method.lower() == 'fortran':
        lamda, a, nexp = NemohCore.old_prony_decomposition.lisc(dimensionless_omega, dimensionless_wavenumber)
        lamda = lamda[:nexp]
        a = a[:nexp]

    else:
        raise ValueError("Unrecognized method name for the Prony decomposition.")

    # Add one more exponential function (actually a constant).
    # It is not clear where it comes from exactly in the theory...
    a = np.concatenate([a, np.array([2])])
    lamda = np.concatenate([lamda, np.array([0.0])])

    return a, lamda


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
        mean square error of the decompostion
    """
    a = np.asarray(a)[:, np.newaxis]
    lamda = np.asarray(lamda)[:, np.newaxis]

    def f(x):
        return np.sum(a * np.exp(lamda*x), axis=0)

    return np.square(f(X) - F).mean()

