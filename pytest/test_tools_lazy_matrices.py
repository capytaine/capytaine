import numpy as np

from capytaine.tools.lazy_matrices import LazyMatrix

RNG = np.random.default_rng()

def test_simple_constructor():
    def constructor(sl):
        return np.arange(sl.start, sl.stop).reshape(-1, 1)
    A = LazyMatrix(constructor, shape=(100, 1), chunk_size=5)
    assert np.all(np.array(A) == np.arange(0, 100).reshape(-1, 1))

def test_mvp():
    def constructor(sl):
        return np.ones((sl.stop-sl.start, 10))
    lazy_ones = LazyMatrix(constructor, shape=(100, 10), chunk_size=5)
    x = RNG.uniform(size=10)
    y1 = lazy_ones @ x
    full_ones = np.ones(lazy_ones.shape)
    y2 = full_ones @ x
    assert np.allclose(y1, y2)
