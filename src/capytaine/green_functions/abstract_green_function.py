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
"""Abstract structure of a class used to compute the Green function"""

from abc import ABC, abstractmethod

import numpy as np


class GreenFunctionEvaluationError(Exception):
    pass


class AbstractGreenFunction(ABC):
    """Abstract method to evaluate the Green function."""

    floating_point_precision: str

    def _get_colocation_points_and_normals(self, mesh1, mesh2, adjoint_double_layer):
        try:
            collocation_points = mesh1.faces_centers
            nb_collocation_points = mesh1.nb_faces
            if not adjoint_double_layer: # Computing the D matrix
                early_dot_product_normals = mesh2.faces_normals
            else: # Computing the K matrix
                early_dot_product_normals = mesh1.faces_normals
            return collocation_points, early_dot_product_normals
        except AttributeError:
            pass

        if isinstance(mesh1, np.ndarray) and mesh1.ndim == 2 and mesh1.shape[1] == 3:
            # This is used when computing potential or velocity at given points in postprocessing
            collocation_points = mesh1
            nb_collocation_points = mesh1.shape[0]
            if not adjoint_double_layer: # Computing the D matrix
                early_dot_product_normals = mesh2.faces_normals
            else: # Computing the K matrix
                early_dot_product_normals = np.zeros((nb_collocation_points, 3))
                # Dummy argument since this method is meant to be used either
                # - to compute potential, then only S is needed and early_dot_product_normals is irrelevant,
                # - to compute velocity, then the adjoint full gradient is needed and early_dot_product is False and this value is unused.
                # TODO: add an only_S argument and return an error here if (early_dot_product and not only_S)
            return collocation_points, early_dot_product_normals

        else:
            raise ValueError(f"Unrecognized first input for {self.__class__.__name__}.evaluate:\n{mesh1}")


    def _init_matrices(self, shape, early_dot_product):
        if self.floating_point_precision == "float32":
            dtype = "complex64"
        elif self.floating_point_precision == "float64":
            dtype = "complex128"
        else:
            raise NotImplementedError(
                    f"Unsupported floating point precision: {self.floating_point_precision}"
                    )

        S = np.zeros(shape, order="F", dtype=dtype)
        K = np.zeros((1 if early_dot_product else 3, shape[0], shape[1]), order="F", dtype=dtype)
        return S, K

    @abstractmethod
    def evaluate(
        self,
        mesh1,
        mesh2,
        free_surface,
        water_depth,
        wavenumber,
        adjoint_double_layer=True,
        early_dot_product=True,
        diagonal_term_in_double_layer=True,
    ):
        pass
