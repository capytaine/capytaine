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
from capytaine.tools.symbolic_multiplication import SymbolicMultiplication, supporting_symbolic_multiplication

def test_definition():
    zero = SymbolicMultiplication("0")
    assert zero.symbol == "0"
    assert isinstance(zero, SymbolicMultiplication)

def test_multiplication():
    zero = SymbolicMultiplication("0")
    b = 2 * zero
    assert b.value == 2 * zero.value

def test_division():
    zero = SymbolicMultiplication("0")
    b = 2 * zero
    assert b / zero == 2

def test_double_division():
    zero = SymbolicMultiplication("0")
    b = 7 * zero * zero
    assert b / (zero * zero) == 7

def test_invert_division():
    zero = SymbolicMultiplication("0")
    b = 9 * zero
    assert zero / b == pytest.approx(1/9)

def test_float():
    zero = SymbolicMultiplication("0")
    b = 4.5 * zero
    assert float(b) == pytest.approx(0.0)

def test_numpy_array():
    zero = SymbolicMultiplication("0")
    b = np.random.rand(10) * zero
    assert (b / zero).shape == (10,)

def test_numpy_array_sum():
    zero = SymbolicMultiplication("0")
    b = np.ones(10) * zero
    assert (np.sum(b) / zero) == 10

def test_numpy_array_sum_with_initialization():
    zero = SymbolicMultiplication("0")
    b = np.ones(10) * zero
    with pytest.raises(TypeError):  # Not implemented (yet?) for these types
        np.sum(b, zero*0)

def test_numpy_matmul():
    zero = SymbolicMultiplication("0")
    b = np.random.rand(10) * zero
    A = np.random.rand(10, 10)
    c = A @ b
    assert (c/zero).shape == (10,)

def test_numpy_einsum():
    zero = SymbolicMultiplication("0")
    b = np.random.rand(10) * zero
    A = np.random.rand(10, 10)
    c = np.einsum('ij,j->i', A, b)
    assert (c/zero).shape == (10,)

def test_numpy_concatenate():
    zero = SymbolicMultiplication("0")
    a = np.ones((10,)) * zero
    b = np.zeros((10,)) * zero
    ab = np.concatenate([a, b], axis=0)
    assert ab.shape == (20,)

def test_numpy_stack():
    zero = SymbolicMultiplication("0")
    a = np.ones((10,)) * zero
    b = np.zeros((10,)) * zero
    ab = np.stack([a, b], axis=0)
    assert ab.shape == (2, 10)

def test_comparison():
    a = SymbolicMultiplication("0", np.arange(10))
    assert np.all(a == 0)

def test_setitem_mask():
    a = SymbolicMultiplication("0", np.arange(10))
    b = SymbolicMultiplication("0", 99)
    a[a.value > 5] = b
    assert np.all(a.value == np.array([0, 1, 2, 3, 4, 5, 99, 99, 99, 99]))

def test_undefined_case():
    assert np.isnan(float(SymbolicMultiplication("0", np.inf)))
    assert np.isnan(float(SymbolicMultiplication("∞", 0.0)))

def test_supporting_symbolic_multiplication():
    zero = SymbolicMultiplication("0")

    @supporting_symbolic_multiplication
    def my_linear_operator(A, x):
        return np.linalg.solve(A, x)

    b = np.random.rand(10) * zero
    A = np.random.rand(10, 10)
    c = my_linear_operator(A, b)
    assert (c/zero).shape == (10,)
