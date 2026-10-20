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


def _field_axis(axis):
    """Check that `axis` (int or tuple of int) is negative: the field axes are at the end."""
    for a in (axis if isinstance(axis, tuple) else (axis,)):
        if a >= 0:
            raise ValueError(f"Use a negative axis instead of {a}: the field axes are at the end, after the batch axes.")
    return axis


class PerturbationExpansion:
    """A quantity x0 + eps*x1 + eps^2*x2 + ... + eps^N*xN, truncated at order N.

    Each coefficient may be a plain scalar or a numpy array (vector, matrix, ...).
    Use `.coefficients[k]` to get the raw coefficient of a specific order.
    Two expansions can only be combined if they share the same truncation order N.
    Plain scalars and arrays used as operands are constants: `x` stands for `x + 0*eps + ...`.

    The axes of the "field" (e.g. `(nb_faces, 3)` for a vector field) are at the end of the coefficients.
    Leading axes are batch axes, such as the wave axes of `DifferenceFrequencyExpansion`.
    Operations act on each coefficient with the usual NumPy broadcasting rules.
    Methods taking an axis (`sum`, `expand_dims`, `squeeze`) take a negative axis (`axis=-1` for the last one),
    and `trace` and `T` act on the last two axes.
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

    def _new(self, coefficients):
        """New expansion of the same class as `self`."""
        return type(self)(coefficients)

    def _map(self, f):
        """Apply `f` to each coefficient (for linear operations)."""
        return self._new([f(c) for c in self.coefficients])

    def _as_expansion(self, other):
        """`other` as an expansion that can be combined with `self`."""
        if isinstance(other, PerturbationExpansion):
            if type(other) is not type(self):
                raise TypeError(f"Cannot combine a {type(self).__name__} with a {type(other).__name__}.")
            _check_same_order(self, other)
            return other
        return self.constant(other, self.order)

    def _aligned(self, other):
        """Coefficients of `self` and `other`, ready to be combined order by order."""
        other = self._as_expansion(other)
        return self.coefficients, other.coefficients

    def _cauchy_product(self, a, b, product):
        """Truncated product of the series of coefficients `a` and `b`."""
        result = []
        for k in range(self.order + 1):
            total = product(a[0], b[k])
            for i in range(1, k + 1):
                total = total + product(a[i], b[k - i])
            result.append(total)
        return result

    def __getitem__(self, index):
        return self._map(lambda c: c[index])

    def __len__(self):
        return len(self.coefficients[0])  # Assumed to be the same for all other coefficients[i]

    def __repr__(self):
        return f"{type(self).__name__}(order={self.order}, coefficients={self.coefficients!r})"

    def __add__(self, other):
        a, b = self._aligned(other)
        return self._new([x + y for x, y in zip(a, b)])

    __radd__ = __add__

    def __neg__(self):
        return self._map(lambda c: -c)

    def __sub__(self, other):
        a, b = self._aligned(other)
        return self._new([x - y for x, y in zip(a, b)])

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        """Elementwise (Hadamard/scalar) truncated Cauchy product."""
        if isinstance(other, (int, float, complex)):
            return self._new([other * c for c in self.coefficients])
        else:
            a, b = self._aligned(other)
            return self._new(self._cauchy_product(a, b, operator.mul))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, float, complex, np.ndarray)):
            return self._new([c / other for c in self.coefficients])
        else:
            return NotImplemented

    def __matmul__(self, other):
        """Matrix-multiplication truncated Cauchy product (for matrix/vector-valued coefficients)."""
        a, b = self._aligned(other)
        return self._new(self._cauchy_product(a, b, np.matmul))

    def __rmatmul__(self, other):
        a, b = self._aligned(other)
        return self._new(self._cauchy_product(b, a, np.matmul))

    def sum(self, axis):
        """Sum over the field axis (or tuple of axes) `axis` (negative)."""
        axis = _field_axis(axis)
        return self._map(lambda c: np.sum(c, axis=axis))

    def expand_dims(self, axis):
        """Insert a field axis of size 1 at position `axis` (negative)."""
        axis = _field_axis(axis)
        return self._map(lambda c: np.expand_dims(c, axis=axis))

    def squeeze(self, axis):
        """Remove the field axis of size 1 at position `axis` (negative)."""
        axis = _field_axis(axis)
        return self._map(lambda c: np.squeeze(c, axis=axis))

    def trace(self):
        """Trace over the last two field axes."""
        return self._map(lambda c: np.trace(c, axis1=-2, axis2=-1))

    def transpose(self):
        """Swap the last two field axes."""
        return self._map(lambda c: np.swapaxes(c, -1, -2))

    @property
    def T(self):
        return self.transpose()

    def _identity_and_perturbation(self, name):
        """Identity matrix I and expansion A such that self = I + A, A having no order 0."""
        if self.order > 2:
            raise NotImplementedError(f"Matrix {name} is only implemented for perturbation of order up to 2")
        I = np.eye(np.shape(self.coefficients[0])[-1])
        A = self - I
        if not np.all(np.abs(A.coefficients[0]) <= 1e-10):
            raise NotImplementedError(f"Matrix {name} is only implemented for perturbation of identity matrix")
        return I, A

    def inv(self):
        """Matrix inverse. Expansion of (I + A)^-1 = I - A + A @ A - ..., up to order 2."""
        I, A = self._identity_and_perturbation("inverse")
        return I - A + A @ A

    def det(self):
        """Matrix determinant. Expansion of det(I + A) = 1 + tr(A) + 1/2 (tr(A)^2 - tr(A @ A)) + ..., up to order 2."""
        _, A = self._identity_and_perturbation("determinant")
        trA = A.trace()
        return 1.0 + trA + 0.5 * (trA * trA - (A @ A).trace())

    def evaluate(self, eps):
        """Numerically sum the series at a concrete eps (mainly for validation/testing)."""
        result = self.coefficients[0]
        power = 1.0
        for k in range(1, self.order + 1):
            power = power * eps
            result = result + power * self.coefficients[k]
        return result


class DifferenceFrequencyExpansion(PerturbationExpansion):
    """Order-2 expansion of a field in a sea made of several wave components at one frequency.

    For a field of shape S, with W wave components (e.g. wave directions):

    - order 0, shape `S`: still water value (geometry, hydrostatic pressure, ...). Must be real.
    - order 1, shape `(W, *S)`: `x1[m]` is the complex amplitude of the response to wave m alone.
    - order 2, shape `(W, W, *S)`: `x2[m, n]` is the difference-frequency QTF: the second-order
      part of the field is `sum_mn x2[m, n] exp(-i (omega_m - omega_n) t)`. Hermitian in (m, n).

    Products keep the difference-frequency part of the product of first-order terms:
    `Q_mn = 1/4 (a_m conj(b_n) + conj(a_n) b_m)`, such that `Q_mm = 1/2 Re(a_m conj(b_m))` is the mean value.
    A complex factor (e.g. i*omega) must be applied to the first-order amplitudes before
    building the expansion, never to an expansion that has an order-2 part.

    The order-0 and order-2 coefficients can be given as scalars (e.g. 0), and are then broadcast.
    When two operands have a different number of field axes, size-1 axes are inserted in front
    of the field axes of the one that has fewer (NumPy broadcasting of the field axes only).
    The wave axes are leading batch axes and the field axes are at the end: index with `x[..., i]`.
    """

    def __init__(self, coefficients):
        x0, x1, x2 = (np.asarray(c) for c in coefficients)
        if np.iscomplexobj(x0):
            if np.max(np.abs(x0.imag)) < 1e-10:
                x0 = x0.real
            else:
                raise ValueError("The order-0 coefficient of a DifferenceFrequencyExpansion must be real. "
                                 "Apply complex factors to the first-order amplitudes before building the expansion.")
        if x1.ndim == 0 or x2.ndim == 1:
            raise ValueError("The order-1 coefficient needs one leading wave axis and the order-2 coefficient two (or none if it is a scalar).")
        nb_waves = x1.shape[0]
        field_shape = np.broadcast_arrays(x0, x1[0], x2 if x2.ndim == 0 else x2[0, 0])[0].shape
        super().__init__([
            np.broadcast_to(x0, field_shape),
            np.broadcast_to(x1, (nb_waves,) + field_shape),
            np.broadcast_to(x2, (nb_waves, nb_waves) + field_shape),
        ])

    @classmethod
    def first_order(cls, amplitudes):
        """First-order quantity: zero order 0 and order 2, and `amplitudes` of shape `(W, *S)` as order 1."""
        return cls([0.0, amplitudes, 0.0])

    @property
    def nb_waves(self):
        return self.coefficients[1].shape[0]

    @property
    def field_ndim(self):
        return self.coefficients[0].ndim

    def constant(self, value, order=2):
        """Embed `value` as the order-0 term, zero elsewhere (stored as broadcast views)."""
        value = np.asarray(value)
        zero = np.zeros(())
        return DifferenceFrequencyExpansion([value, np.broadcast_to(zero, (self.nb_waves,) + value.shape), zero])

    def _as_expansion(self, other):
        other = super()._as_expansion(other)
        if other.nb_waves != self.nb_waves:
            raise ValueError(f"Cannot combine expansions with {self.nb_waves} and {other.nb_waves} waves.")
        return other

    def _with_field_ndim(self, field_ndim):
        """Coefficients with size-1 axes inserted after their k wave axes, up to `field_ndim` field axes."""
        missing = (1,) * (field_ndim - self.field_ndim)
        return [np.reshape(c, c.shape[:k] + missing + c.shape[k:]) for k, c in enumerate(self.coefficients)]

    def _aligned(self, other):
        other = self._as_expansion(other)
        field_ndim = max(self.field_ndim, other.field_ndim)
        return self._with_field_ndim(field_ndim), other._with_field_ndim(field_ndim)

    def _cauchy_product(self, a, b, product):
        return [
            product(a[0], b[0]),
            product(a[0], b[1]) + product(a[1], b[0]),
            product(a[0], b[2]) + self._cross(a[1], b[1], product) + product(a[2], b[0]),
        ]

    @staticmethod
    def _cross(a1, b1, product):
        """Q_mn = 1/4 (a_m conj(b_n) + conj(a_n) b_m), for all (m, n) at once."""
        a_m, b_m = np.expand_dims(a1, axis=1), np.expand_dims(b1, axis=1)  # (W, 1, *S)
        a_n, b_n = np.expand_dims(a1, axis=0), np.expand_dims(b1, axis=0)  # (1, W, *S)
        return (product(a_m, np.conj(b_n)) + product(np.conj(a_n), b_m)) / 4

    def _check_matrix_operand(self, other):
        # Vectors would be broadcast as matrices after the alignment of the field axes.
        other_ndim = other.field_ndim if isinstance(other, PerturbationExpansion) else np.ndim(other)
        if min(self.field_ndim, other_ndim) < 2:
            raise ValueError("Matrix product of a DifferenceFrequencyExpansion needs operands with at least 2 field axes: "
                             "use `expand_dims(-1)` for a vector, or `(a * b).sum(-1)` for a dot product.")

    def __matmul__(self, other):
        self._check_matrix_operand(other)
        return super().__matmul__(other)

    def __rmatmul__(self, other):
        self._check_matrix_operand(other)
        return super().__rmatmul__(other)

    def __getitem__(self, index):
        if not (isinstance(index, tuple) and len(index) > 0 and index[0] is Ellipsis):
            raise IndexError("The field axes are at the end, after the wave axes: index with `x[..., i]`.")
        return super().__getitem__(index)

    def evaluate(self, eps):
        raise NotImplementedError("The coefficients of a DifferenceFrequencyExpansion are amplitudes at different frequencies "
                                  "(omega_m at order 1, omega_m - omega_n at order 2): their sum at a given eps has no meaning.")
