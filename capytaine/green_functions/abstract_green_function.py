"""Abstract structure of a class used to compute the Green function"""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC, abstractmethod

import numpy as np

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes

class AbstractGreenFunction(ABC):
    """Abstract method to evaluate the Green function."""

    def _get_colocation_points_and_normals(self, mesh1, mesh2, adjoint_double_layer):
        if isinstance(mesh1, Mesh) or isinstance(mesh1, CollectionOfMeshes):
            collocation_points = mesh1.faces_centers
            nb_collocation_points = mesh1.nb_faces
            if not adjoint_double_layer: # Computing the D matrix
                early_dot_product_normals = mesh2.faces_normals
            else: # Computing the K matrix
                early_dot_product_normals = mesh1.faces_normals

        elif isinstance(mesh1, np.ndarray) and mesh1.ndim == 2 and mesh1.shape[1] == 3:
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

        else:
            raise ValueError(f"Unrecognized first input for {self.__class__.__name__}.evaluate:\n{mesh1}")

        return collocation_points, early_dot_product_normals

    def _init_matrices(self, shape, dtype, early_dot_product):
        S = np.empty(shape, order="F", dtype=dtype)
        K = np.empty((shape[0], shape[1], 1 if early_dot_product else 3), order="F", dtype=dtype)
        return S, K

    @abstractmethod
    def evaluate(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=True):
        pass
