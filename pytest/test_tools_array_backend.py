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
import pytest

import numpy as np
xps = pytest.importorskip("array_api_strict")  # Not available on Python 3.8

from capytaine.tools.array_backend import is_array, complex_dtype, to_backend_of, to_numpy
from capytaine.tools.symbolic_multiplication import SymbolicMultiplication


def test_is_array():
    assert is_array(np.zeros(3))
    assert is_array(xps.zeros(3))
    assert not is_array([1.0, 2.0])
    assert not is_array(1.0)
    assert not is_array(SymbolicMultiplication("0", np.zeros(3)))


@pytest.mark.parametrize("xp", [np, xps])
def test_complex_dtype_preserves_precision(xp):
    assert complex_dtype(xp, xp.float32) == xp.complex64
    assert complex_dtype(xp, xp.complex64) == xp.complex64
    assert complex_dtype(xp, xp.float64) == xp.complex128
    assert complex_dtype(xp, xp.complex128) == xp.complex128


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("complex_", [False, True])
def test_to_backend_of_numpy_to_strict(dtype, complex_):
    like = xps.zeros(3, dtype=getattr(xps, dtype))
    x = to_backend_of(np.array([1.0, 2.0, 3.0]), like, complex_=complex_)
    assert isinstance(x, type(like))
    expected = {("float32", False): xps.float32, ("float32", True): xps.complex64,
                ("float64", False): xps.float64, ("float64", True): xps.complex128}
    assert x.dtype == expected[(dtype, complex_)]
    np.testing.assert_allclose(to_numpy(x), [1.0, 2.0, 3.0])


def test_to_backend_of_strict_to_numpy():
    x = to_backend_of(xps.asarray([1.0, 2.0]), np.zeros(2, dtype=np.float32), complex_=True)
    assert isinstance(x, np.ndarray)
    assert x.dtype == np.complex64


@pytest.mark.parametrize("make", [
    lambda: np.array([1.0, 2.0]),
    lambda: xps.asarray([1.0, 2.0]),
    lambda: xps.asarray([1.0 + 2.0j, 2.0 - 1.0j]),
])
def test_to_numpy_round_trip(make):
    x = make()
    y = to_numpy(x)
    assert isinstance(y, np.ndarray)
    np.testing.assert_array_equal(y, np.asarray(x))


def test_to_numpy_of_none_and_numpy_array():
    assert to_numpy(None) is None
    x = np.zeros(3)
    assert to_numpy(x) is x


def test_to_numpy_of_symbolic_multiplication():
    s = SymbolicMultiplication("0", xps.asarray([1.0, 2.0], dtype=xps.complex64))
    back = to_numpy(s)
    assert isinstance(back, SymbolicMultiplication)
    assert back.symbol == "0"
    assert isinstance(back.value, np.ndarray)
    assert back.value.dtype == np.complex64
    np.testing.assert_allclose(back.value, [1.0, 2.0])
