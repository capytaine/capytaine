"""Prony decomposition: tool to approximate a function as a sum of exponentials.
Used in particular in the finite depth Green function.
"""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

import numpy as np
from numpy.polynomial import polynomial
from scipy.optimize import curve_fit
from scipy.linalg import toeplitz

LOG = logging.getLogger(__name__)
RNG = np.random.default_rng()


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
    """Mean square error of the exponential decomposition defined by the
    coefficients a and lamda with respect to the reference values in F.

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


class PronyDecompositionFailure(Exception):
    pass


def find_best_exponential_decomposition(f, x_min, x_max, n_exp_range, *, tol=1e-4, noise_on_domain_points_std=0.01):
    """Tries to construct an exponential decompositoin of the function f on the
    domain [x_min, x_max] by testing the number of exponentials in n_exp_range.

    Parameters
    ----------
    f: callable
        The function ℝ→ℝ to be approximated.
        Should support vectorized calls (that is passing a vector of inputs
        and get the vector of corresponding outputs)
    x_min, x_max: floats
        The bounds of the domain of input in which f should be approximated
    n_exp_range: iterable of ints
        The decomposition sizes that will be tested
    tol: float, optional
        The target mean square error.
    noise_on_domain_points_std: float, optional
        Introduces some random variability on the points where the function is evaluated.
        Set this parameter to zero to disable randomness.

    """
    # Try different range of evaluation points to construct the decomposition.
    for n_exp in n_exp_range:

        # f might be ill-defined at some single specific values
        # (for the use-case of delhommeau.py, it is when x = kh exactly).
        # Thus we slightly randomize the range of evaluation points for the Prony decomposition.
        # This way, if one of the evaluation points hits the singular point, it will most likely not hit it again at the next iteration.
        x_max_iter = (1 + noise_on_domain_points_std*RNG.uniform())*x_max

        try:
            # The coefficients are computed on a resolution of 4*n_exp+1 ...
            X = np.linspace(x_min, x_max_iter, 4*n_exp+1)
            a, lamda = exponential_decomposition(X, f(X), n_exp)

            # ... and they are evaluated on a finer discretization.
            X = np.linspace(x_min, x_max_iter, 8*n_exp+1)
            if error_exponential_decomposition(X, f(X), a, lamda) < tol:
                return a, lamda
        except Exception:
            # If something bad happened while computing the decomposition, try
            # the next one.
            continue

    raise PronyDecompositionFailure(
            "No suitable Prony decomposition has been found in "
            f"[{x_min}, {x_max}] for tol={tol} "
            f"using a number of terms in {n_exp_range}."
            )
