#!/usr/bin/env python
# coding: utf-8

from collections import OrderedDict
from functools import wraps

def delete_first_lru_cache(maxsize=1):
    """Behaves like functools.lru_cache(), but the oldest data in the cache is
    deleted *before* computing a new one, in order to limit RAM usage when
    stored objects are big."""

    def decorator(f):
        cache = OrderedDict()

        @wraps(f)
        def decorated_f(*args, **kwargs):
            # /!\ cache only args

            if args in cache:
                # Get item in cache
                return cache[args]

            if len(cache) + 1 > maxsize:
                # Drop oldest item in cache.
                cache.popitem(last=False)

            # Compute and store
            result = f(*args, **kwargs)
            cache[args] = result

            return result

        return decorated_f

    return decorator


# if __name__ == "__main__":
#     import numpy as np
#     from functools import lru_cache

#     # @lru_cache(maxsize=1)
#     @delete_first_lru_cache(maxsize=1)
#     def f(i):
#         print(i)
#         np.random.seed(i)
#         a = np.random.random((20_000, 20_000))
#         print(a.nbytes / 2**30, "Go")
#         return a

#     print(f(1)[0, 0])
#     print(f(1)[0, 0])
#     print(f(2)[0, 0])
#     print(f(1)[0, 0])
#     print(f(3)[0, 0])
