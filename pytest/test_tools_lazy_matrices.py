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
