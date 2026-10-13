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
"""Truncated perturbation series ("higher-order dual numbers")."""

import operator

import numpy as np


def _check_same_order(a, b):
    if a.order != b.order:
        raise ValueError(
            f"Cannot combine PerturbationExpansion of order {a.order} "
            f"with one of order {b.order}: orders must match exactly."
        )


def _as_expansion(value, order):
    if isinstance(value, PerturbationExpansion):
        return value
    return PerturbationExpansion.constant(value, order)


def _cauchy_product(a_coefficients, b_coefficients, order, product):
    result = []
    for k in range(order + 1):
        total = product(a_coefficients[0], b_coefficients[k])
        for i in range(1, k + 1):
            total = total + product(a_coefficients[i], b_coefficients[k - i])
        result.append(total)
    return result


class PerturbationExpansion:
    """A quantity x0 + eps*x1 + eps^2*x2 + ... + eps^N*xN, truncated at order N.

    Each coefficient may be a plain scalar or a numpy array (vector, matrix, ...).
    Use `.coefficients[k]` to get the raw coefficient of a specific order.
    Two expansions can only be combined if they share the same truncation order N.
    """

    # Ensure numpy defers to our __radd__/__rmul__/etc. when a bare ndarray is on
    # the left of an operator with a PerturbationExpansion on the right.
    __array_priority__ = 1000

    def __init__(self, coefficients):
        coefficients = tuple(coefficients)
        if len(coefficients) == 0:
            raise ValueError("A PerturbationExpansion needs at least an order-0 coefficient.")
        self.coefficients = coefficients
        self.order = len(coefficients) - 1

    @classmethod
    def constant(cls, value, order):
        """Embed `value` as the order-0 term of an order-`order` expansion, zero elsewhere."""
        if order < 0:
            raise ValueError("order must be >= 0")
        zero = value * 0
        return cls([value] + [zero for _ in range(order)])

    def __getitem__(self, index):
        return PerturbationExpansion([c[index] for c in self.coefficients])

    def __len__(self):
        return len(self.coefficients[0])  # Assumed to be the same for all other coefficients[i]

    def __repr__(self):
        return f"PerturbationExpansion(order={self.order}, coefficients={self.coefficients!r})"

    def __add__(self, other):
        other = _as_expansion(other, self.order)
        _check_same_order(self, other)
        return PerturbationExpansion(
            [a + b for a, b in zip(self.coefficients, other.coefficients)]
        )

    __radd__ = __add__

    def __neg__(self):
        return PerturbationExpansion([-c for c in self.coefficients])

    def __sub__(self, other):
        other = _as_expansion(other, self.order)
        _check_same_order(self, other)
        return PerturbationExpansion(
            [a - b for a, b in zip(self.coefficients, other.coefficients)]
        )

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        """Elementwise (Hadamard/scalar) truncated Cauchy product."""
        other = _as_expansion(other, self.order)
        _check_same_order(self, other)
        coefficients = _cauchy_product(self.coefficients, other.coefficients, self.order, operator.mul)
        return PerturbationExpansion(coefficients)

    __rmul__ = __mul__

    def __matmul__(self, other):
        """Matrix-multiplication truncated Cauchy product (for matrix/vector-valued coefficients)."""
        other = _as_expansion(other, self.order)
        _check_same_order(self, other)
        coefficients = _cauchy_product(self.coefficients, other.coefficients, self.order, np.matmul)
        return PerturbationExpansion(coefficients)

    def __rmatmul__(self, other):
        other = _as_expansion(other, self.order)
        _check_same_order(self, other)
        coefficients = _cauchy_product(other.coefficients, self.coefficients, self.order, np.matmul)
        return PerturbationExpansion(coefficients)

    def transpose(self):
        return PerturbationExpansion([np.swapaxes(c, -1, -2) for c in self.coefficients])

    @property
    def T(self):
        return self.transpose()

    def inv(self):
        """Matrix inverse. Expansion of (I + εA1 + ε²A2)^-1)"""
        n = self.coefficients[0].shape[-1]
        if not np.allclose(self.coefficients[0], np.eye(n)):
            raise NotImplementedError("Matrix determinant is only implemented for perturbation of identity matrix")

        if self.order == 0:
            return PerturbationExpansion([np.eye(n)])
        elif self.order == 1:
            A1 = self.coefficients[1]
            return PerturbationExpansion([np.eye(n), -A1])
        elif self.order == 2:
            A1 = self.coefficients[1]
            A2 = self.coefficients[2]
            return PerturbationExpansion([np.eye(n), -A1, -A2 + A1**2])
        else:
            raise NotImplementedError("Matrix determinant is only implemented for perturbation of order up to 2")

    def det(self):
        """Matrix determinant. Expansion of det(I + εA1 + ε²A2)"""
        if not np.allclose(self.coefficients[0], np.eye(self.coefficients[0].shape[-1])):
            raise NotImplementedError("Matrix determinant is only implemented for perturbation of identity matrix")

        if self.order == 0:
            return PerturbationExpansion([1.0])
        elif self.order == 1:
            A1 = self.coefficients[1]
            return PerturbationExpansion([1.0, np.trace(A1)])
        elif self.order == 2:
            A1 = self.coefficients[1]
            A2 = self.coefficients[2]
            return PerturbationExpansion([1.0, np.trace(A1), np.trace(A2) + 1/2 * (np.trace(A1)**2 - np.trace(A1**2))])
        else:
            raise NotImplementedError("Matrix determinant is only implemented for perturbation of order up to 2")

    def evaluate(self, eps):
        """Numerically sum the series at a concrete eps (mainly for validation/testing)."""
        result = self.coefficients[0]
        power = 1.0
        for k in range(1, self.order + 1):
            power = power * eps
            result = result + power * self.coefficients[k]
        return result
