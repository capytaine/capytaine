from functools import lru_cache

import pytest
import numpy as np

from capytaine.tools.lru_cache import lru_cache_with_strict_maxsize


def get_random_id(n):
    # An impure function that shouldn't be memoized, but allows to test the cache.
    return np.random.default_rng().integers(n)

N = 100_000_000

def test_cache():
    f = lru_cache_with_strict_maxsize()(get_random_id)
    a = f(N)
    b = f(N)
    assert a == b  # randint was not run the second time, the cached value was used


def test_cache_maxsize():
    # 1/N chance of failing
    f = lru_cache_with_strict_maxsize(maxsize=1)(get_random_id)
    a = f(N)
    f(10)
    b = f(N)
    assert a != b  # the cached value of the first call had been deleted


def test_cache_kwargs():
    f = lru_cache_with_strict_maxsize()(get_random_id)
    a = f(n=N)
    b = f(n=N)
    assert a == b  # the cached value of the first call had been deleted


def test_delete_first():
    """This test is a convoluted way to test the main difference between the
    built-in lru_cache and the lru_cache_with_strict_maxsize decorator defined
    by Capytaine. The difference is that the latter make sure that objects in
    the cache are deleted before computing and caching a new result. For
    Capytaine, it is meant to limit the RAM usage. Here, it is tested by making
    sure that there is never two instances of a Singleton class at the same
    time."""

    class Singleton:
        nb_instances = [0]
        # Use a one-element list as a mutable integer shared by all instances.

        def __init__(self):
            if self.nb_instances[0] != 0:
                raise ValueError("There can be only one Singleton at the same time!")
            self.nb_instances[0] += 1

        def __del__(self):
            # When the singleton is deleted, update the counter
            self.nb_instances[0] -= 1

    assert Singleton.nb_instances[0] == 0
    a = Singleton()
    assert Singleton.nb_instances[0] == 1
    del a
    assert Singleton.nb_instances[0] == 0

    ## NO CACHE
    def new_instance(a):
        return Singleton()
    new_instance(1)
    assert Singleton.nb_instances[0] == 0  # The Singleton created above is immediately garbage-collected.

    ## CAPYTAINE'S CACHE
    @lru_cache_with_strict_maxsize(maxsize=1)
    def new_instance(a):
        return Singleton()

    new_instance(1)
    assert Singleton.nb_instances[0] == 1  # It is not garbage collected only because it is still in the cache.
    new_instance(2)  # When we create a new instance, the old one is removed from the cache before creating the new one.

    del new_instance  # Delete the cache before doing more tests.

    # STDLIB CACHE
    @lru_cache(maxsize=1)
    def new_instance(a):
        return Singleton()

    new_instance(1)
    with pytest.raises(ValueError):
        new_instance(2)  # lru_cache tries to create the new singleton, before deleting the old one from the cache.


# def benchmark_ram_usage():
#     # Would need a way to monitor the RAM usage to automate this test
#     # The RAM usage should never be more than the ~3GB allocated by a single call
#
#     # @delete_first_lru_cache(maxsize=1)
#     @lru_cache(maxsize=1)
#     def f(i):
#         np.random.seed(i)
#         a = np.random.random((20_000, 20_000))
#         print(f"Allocated {a.nbytes / 2**30:.2f} GB")
#         return a
#
#     for i in (1, 1, 2, 1, 3):
#         print(f"Calling f({i})")
#         f(i)
