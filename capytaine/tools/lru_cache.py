# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
"""Tools for memoization of functions."""
from collections import OrderedDict
from functools import wraps

import logging

LOG = logging.getLogger(__name__)


def lru_cache_with_strict_maxsize(maxsize=1):
    """Behaves mostly like functools.lru_cache(), but the oldest data in the cache is
    deleted *before* computing a new one, in order to *never* have more that
    `maxsize` items in memory.
    This is useful to limit RAM usage when stored objects are big, like the interaction
    matrices of Capytaine."""

    def decorator(f):
        cache = OrderedDict()

        @wraps(f)
        def decorated_f(*args, **kwargs):
            hashable_kwargs = tuple((k, v) for (k, v) in kwargs.items())
            # Might miss a cache hit if the order of kwargs is changed.
            # But at least unlike a previous version, should not return a wrong value.

            if (args, hashable_kwargs) in cache:
                # Get item in cache
                LOG.debug("Get cached version of %s(%s, %s)", f.__name__, args, hashable_kwargs)
                return cache[(args, hashable_kwargs)]

            if len(cache) + 1 > maxsize:
                # Drop oldest item in cache.
                cache.popitem(last=False)

            # Compute and store
            LOG.debug("Computing %s(%s, %s)", f.__name__, args, hashable_kwargs)
            result = f(*args, **kwargs)
            cache[(args, hashable_kwargs)] = result

            return result

        return decorated_f

    return decorator


delete_first_lru_cache = lru_cache_with_strict_maxsize  # For backward compatibility...
