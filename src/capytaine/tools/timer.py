# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""A simple timer class used to measure the time spent in various parts of the BEM solver."""

from functools import wraps
import time
import contextlib

import pandas as pd

class Timer:
    """A timer class that can be used as context manager or as decorator using `wraps_function` method.
    Several timing measurement can be nested.

    Attributes
    ----------
    timings: List[Dict]]
        List of records of each timing measurement.
        The record is a dict with a 'timing' key and any number of other metadata keys.
    default_tags: Optional[Dict]
        Tags added to all the timing measurements.
    _start_times: List[float]
        Start times of the ongoing timing measurements.

    Example
    -------
    ::

        from time import sleep  # For testing

        timer = Timer()

        with timer(tag="run 1"):
            sleep(1.0)

        print(timer.total)  # 1.0...

        @timer.wraps_function(tag="run function")
        def my_function():
            sleep(0.5)

        my_function()
        print(timer.total)  # 1.5...
        my_function()
        print(timer.total)  # 2.0...

        with timer(tag="outer"):
            sleep(0.3)
            with timer(tag="inner"):
                sleep(0.3)
            sleep(0.3)

    """
    def __init__(self, *, default_tags=None):
        self.timings = []
        if default_tags is None:
            default_tags = {}
        self.default_tags = default_tags
        self._start_times = []  # Starting time of ongoing sub-timers

    def __repr__(self):
        return f"Timer({self.timings})"

    def add_data_from_other_timer(self, other, **supplementary_tags):
        self.timings.extend([{**t, **supplementary_tags} for t in other.timings])

    @contextlib.contextmanager
    def __call__(self, **tags):
        self._start_times.append(time.perf_counter())
        try:
            yield
        finally:
            timing = time.perf_counter() - self._start_times.pop()
            self.timings.append({'timing': timing, **tags, **self.default_tags})

    def wraps_function(self, **tags):
        def wrapper(f):
            @wraps(f)
            def wrapped_f(*args, **kwargs):
                with self(**tags):
                    out = f(*args, **kwargs)
                return out
            return wrapped_f
        return wrapper

    def as_dataframe(self):
        if len(self.timings) == 0:
            return pd.DataFrame([{'timing': 0.0}])
        else:
            return pd.DataFrame(self.timings)

    @property
    def total(self):
        return self.as_dataframe()['timing'].sum()
