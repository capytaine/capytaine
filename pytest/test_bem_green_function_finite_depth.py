import pytest

import numpy as np
import capytaine as cpt

def test_prony_decomposition():
    gf = cpt.Delhommeau(finite_depth_prony_decomposition_method="python")
    decomp_default = gf.find_best_exponential_decomposition(1.0)
    decomp_f = gf.find_best_exponential_decomposition(1.0, method="fortran")
    decomp_p = gf.find_best_exponential_decomposition(1.0, method="python")
    assert np.allclose(decomp_default[0], decomp_p[0])
    assert np.allclose(decomp_p[0], decomp_f[0], rtol=0.2)
    assert np.allclose(decomp_p[1], decomp_f[1], rtol=0.2)
