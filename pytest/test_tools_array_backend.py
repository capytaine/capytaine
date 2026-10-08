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

try:
    import array_api_strict as xps
except ImportError:  # Not available on Python 3.8
    xps = None

from capytaine.tools.array_backend import (
    IterableTensor, array_namespace, complex_dtype, ending_dimensions_at_the_beginning, is_array,
    is_numpy_namespace, leading_dimensions_at_the_end, split, to_backend_of, to_numpy,
)
from capytaine.tools.symbolic_multiplication import SymbolicMultiplication

requires_strict = pytest.mark.skipif(xps is None, reason="array_api_strict is not available")

RNG = np.random.default_rng(seed=0)


@requires_strict
def test_is_array():
    assert is_array(np.zeros(3))
    assert is_array(xps.zeros(3))
    assert not is_array([1.0, 2.0])
    assert not is_array(1.0)
    assert not is_array(SymbolicMultiplication("0", np.zeros(3)))


@pytest.mark.parametrize("xp", [np, pytest.param(xps, marks=requires_strict)])
def test_complex_dtype_preserves_precision(xp):
    assert complex_dtype(xp, xp.float32) == xp.complex64
    assert complex_dtype(xp, xp.complex64) == xp.complex64
    assert complex_dtype(xp, xp.float64) == xp.complex128
    assert complex_dtype(xp, xp.complex128) == xp.complex128


@requires_strict
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


@requires_strict
def test_to_backend_of_strict_to_numpy():
    x = to_backend_of(xps.asarray([1.0, 2.0]), np.zeros(2, dtype=np.float32), complex_=True)
    assert isinstance(x, np.ndarray)
    assert x.dtype == np.complex64


@pytest.mark.parametrize("make", [
    lambda: np.array([1.0, 2.0]),
    pytest.param(lambda: xps.asarray([1.0, 2.0]), marks=requires_strict),
    pytest.param(lambda: xps.asarray([1.0 + 2.0j, 2.0 - 1.0j]), marks=requires_strict),
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


@requires_strict
def test_to_numpy_of_symbolic_multiplication():
    s = SymbolicMultiplication("0", xps.asarray([1.0, 2.0], dtype=xps.complex64))
    back = to_numpy(s)
    assert isinstance(back, SymbolicMultiplication)
    assert back.symbol == "0"
    assert isinstance(back.value, np.ndarray)
    assert back.value.dtype == np.complex64
    np.testing.assert_allclose(back.value, [1.0, 2.0])


@requires_strict
def test_is_numpy_namespace():
    assert is_numpy_namespace(array_namespace(np.zeros(3)))
    assert not is_numpy_namespace(array_namespace(xps.zeros(3)))


def test_permute_dims():
    a = RNG.normal(size=(1, 2, 3, 4, 5))
    assert leading_dimensions_at_the_end(a).shape == (3, 4, 5, 1, 2)
    assert ending_dimensions_at_the_beginning(a).shape == (4, 5, 1, 2, 3)
    assert np.allclose(ending_dimensions_at_the_beginning(leading_dimensions_at_the_end(a)), a)
    assert np.allclose(leading_dimensions_at_the_end(ending_dimensions_at_the_beginning(a)), a)


@requires_strict
def test_permute_dims_with_another_array_library():
    a = RNG.normal(size=(2, 3, 4))
    np.testing.assert_array_equal(np.asarray(leading_dimensions_at_the_end(xps.asarray(a))), leading_dimensions_at_the_end(a))
    np.testing.assert_array_equal(np.asarray(ending_dimensions_at_the_beginning(xps.asarray(a))), ending_dimensions_at_the_beginning(a))


def test_split():
    x = np.arange(12.0)
    parts = split(x, 4)
    assert len(parts) == 4
    assert all(np.array_equal(p, ref) for p, ref in zip(parts, np.split(x, 4)))


@requires_strict
def test_split_with_another_array_library():
    parts = split(xps.arange(12.0), 3)
    assert [p.shape for p in parts] == [(4,), (4,), (4,)]
    np.testing.assert_array_equal(np.asarray(parts[1]), [4.0, 5.0, 6.0, 7.0])
    # Splitting along the first axis of a 2D array
    parts = split(xps.reshape(xps.arange(12.0), (6, 2)), 3)
    assert [p.shape for p in parts] == [(2, 2), (2, 2), (2, 2)]


def test_split_2d():
    x = np.arange(12.0).reshape((6, 2))
    assert all(np.array_equal(p, ref) for p, ref in zip(split(x, 3), np.split(x, 3)))


def test_iterable_tensor_from_a_list():
    items = [np.full((2, 2), i) for i in range(4)]
    t = IterableTensor(items)
    assert len(t) == 4
    assert all(np.array_equal(a, b) for a, b in zip(t, items))
    a, b, c, d = t
    assert np.array_equal(c, items[2])
    assert np.array_equal(t[-1], items[3])
    assert len(t[1:3]) == 2
    assert t.as_array(np).shape == (4, 2, 2)


def test_iterable_tensor_from_an_array_is_not_copied():
    array = RNG.normal(size=(5, 2, 2))
    t = IterableTensor(array)
    assert len(t) == 5
    assert np.shares_memory(t[2], array)
    assert np.shares_memory(t[1:][0], array)
    assert np.array_equal(t[-1], array[-1])
    assert t.as_array(np) is array
    assert len(list(t)) == 5


@requires_strict
def test_iterable_tensor_from_an_array_of_another_library():
    array = xps.asarray(RNG.normal(size=(3, 2, 2)))
    t = IterableTensor(array)
    assert len(t) == 3
    a, b, c = t  # Not possible with the array itself
    assert b.shape == (2, 2)
    assert [x.shape for x in t[1:]] == [(2, 2), (2, 2)]
    assert t.as_array(xps) is array
