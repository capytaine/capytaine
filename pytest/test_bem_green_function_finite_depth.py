import pytest

import numpy as np
import capytaine as cpt

from capytaine.tools.prony_decomposition import exponential_decomposition, find_best_exponential_decomposition

def test_prony_decomposition():
    x = np.linspace(0.0, 1.0, 100)
    y = 2*np.exp(-x) - 4*np.exp(-2*x)
    a, lamda = exponential_decomposition(x, y, 2)
    assert ((np.allclose(a, [-4.0, 2.0]) and np.allclose(lamda, [-2.0, -1.0]))
            or (np.allclose(a, [2.0, -4.0]) and np.allclose(lamda, [-1.0, -2.0])))


def test_find_best_prony_decomposition():
    rng = np.random.default_rng(seed=987654321)
    # May fail for other RNGs, for instance if two values from lamda_refs are very close from each other
    n_ref = 6
    a_ref = rng.uniform(-10.0, 10.0, size=n_ref)
    lamda_ref = -rng.chisquare(1.0, size=n_ref)
    def f(x):
        return sum(a_ref[i]*np.exp(lamda_ref[i]*x) for i in range(n_ref))
    a, lamda = find_best_exponential_decomposition(f, 0.0, 100.0, range(2, 15), tol=1e-4)
    assert np.allclose(np.sort(a), np.sort(a_ref))
    assert np.allclose(np.sort(lamda_ref), np.sort(lamda_ref))


def test_python_and_fortran_prony_decomposition_for_green_function():
    gf = cpt.Delhommeau(finite_depth_prony_decomposition_method="python")
    decomp_default = gf.find_best_exponential_decomposition(1.0)
    decomp_f = gf.find_best_exponential_decomposition(1.0, method="fortran")
    decomp_p = gf.find_best_exponential_decomposition(1.0, method="python")
    assert np.allclose(decomp_default[0], decomp_p[0])
    assert np.allclose(decomp_p[0], decomp_f[0], rtol=0.2)
    assert np.allclose(decomp_p[1], decomp_f[1], rtol=0.2)
