"""Abstract structure of a class used to compute the Green function"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

from abc import ABC, abstractmethod

class AbstractGreenFunction(ABC):
    """Abstract method to evaluate the Green function."""

    @abstractmethod
    def evaluate(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
        pass
