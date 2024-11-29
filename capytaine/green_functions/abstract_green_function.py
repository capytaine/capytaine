"""Abstract structure of a class used to compute the Green function"""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC, abstractmethod

class AbstractGreenFunction(ABC):
    """Abstract method to evaluate the Green function."""

    @abstractmethod
    def evaluate(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=True):
        pass
