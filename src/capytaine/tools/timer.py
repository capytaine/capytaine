"""A simple timer class used to measure the time spent in various parts of the BEM solver."""

from functools import wraps
import time

class Timer:
    """A simple timer class that can be used as context manager or as decorator using `wraps_function` method

    Example
    -------
    ::

        timer = Timer()
        with timer:
            sleep(1.0)

        print(timer.total)  # 1.0...

        @timer.wraps_function
        def my_function():
            sleep(0.5)

        my_function()
        print(timer.total)  # 1.5...
        my_function()
        print(timer.total)  # 2.0...

        print(timer.timings)  # [1.0, 0.5, 0.5]
    """

    def __init__(self, timings=None):
        if timings is None:
            self.timings = []
        else:
            self.timings = timings

    def __repr__(self):
        return f"Timer({self.timings})"

    @property
    def nb_timings(self):
        return len(self.timings)

    @property
    def total(self):
        return sum(self.timings)

    @property
    def mean(self):
        if self.nb_timings == 0:
            return float('nan')
        else:
            return self.total/self.nb_timings

    def __enter__(self):
        self.start_time = time.perf_counter()

    def __exit__(self, *exc):
        self.timings.append(time.perf_counter() - self.start_time)

    def wraps_function(self, f):
        @wraps(f)
        def wrapped_f(*args, **kwargs):
            with self:
                return f(*args, **kwargs)
        return wrapped_f
