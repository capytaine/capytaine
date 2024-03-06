import numpy as np
from functools import total_ordering, wraps

# @total_ordering
class SymbolicMultiplication:
    def __init__(self, symbol, value=1.0):
        self.symbol = symbol
        self.value = value

    def __format__(self, format_spec):
        return f"{self.symbol}×{self.value.__format__(format_spec)}"

    __array_priority__ = 1.0

    def __array_function__(self, func, types, *args, **kwargs):
        if func in {np.real, np.imag, np.sum}:
            return SymbolicMultiplication(self.symbol, func(self.value))
        else:
            return NotImplemented

    def __str__(self):
        return f"{self.symbol}×{self.value}"

    def __repr__(self):
        return f"SymbolicMultiplication(\"{self.symbol}\", {repr(self.value)})"

    def __mul__(self, x):
        return SymbolicMultiplication(self.symbol, self.value * x)

    def __rmul__(self, x):
        return SymbolicMultiplication(self.symbol, x * self.value)

    def __pow__(self, n):
        if n == 2:
            return self * self
        else:
            raise NotImplementedError

    def __truediv__(self, x):
        if hasattr(x, 'symbol') and self.symbol == x.symbol:
            return self.value / x.value
        else:
            return SymbolicMultiplication(self.symbol, self.value / x)

    def __rtruediv__(self, x):
        if hasattr(x, 'symbol') and self.symbol == x.symbol:
            return x.value / self.value
        elif self.symbol == "0":
            return SymbolicMultiplication("∞", x/self.value)
        elif self.symbol == "∞":
            return SymbolicMultiplication("0", x/self.value)
        else:
            raise NotImplementedError

    def __matmul__(self, x):
        return SymbolicMultiplication(self.symbol, self.value @ x)

    def __rmatmul__(self, x):
        return SymbolicMultiplication(self.symbol, x @ self.value)

    def __eq__(self, x):
        return float(self) == x

    def __lt__(self, x):
        return float(self) < x

    def __hash__(self):
        return hash((self.symbol, self.value))

    def __float__(self):
        if self.symbol == "0":
            return 0.0
        elif self.symbol == "∞":
            return np.inf
        else:
            raise NotImplementedError

    def reshape(self, *args):
        return SymbolicMultiplication(self.symbol, self.value.reshape(*args))


def supporting_symbolic_multiplication(f):
    @wraps(f)
    def wrapped_f(a, x):
        if hasattr(x, 'symbol'):
            return SymbolicMultiplication(x.symbol, f(a, x.value))
        else:
            return f(a, x)
    return wrapped_f
