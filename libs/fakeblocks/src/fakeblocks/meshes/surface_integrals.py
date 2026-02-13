"""Tools for surface integrals and hydrostatics."""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC
import numpy as np

class SurfaceIntegralsMixin(ABC):
    # These methods need to be defined for both Mesh and CollectionOfMeshes with the exact same definitions.
    # To avoid redunduncy, they are defined here in a mixin inherited by both Mesh and CollectionOfMeshes.

    def surface_integral(self, data, **kwargs):
        """Returns integral of given data along wet surface area."""
        return np.sum(data * self.faces_areas, **kwargs)

    @property
    def wet_surface_area(self):
        """Returns wet surface area."""
        return self.immersed_part().surface_integral(1)

    @property
    def volumes(self):
        """Returns volumes using x, y, z components of the mesh."""
        norm_coord = self.faces_normals * self.faces_centers
        return self.surface_integral(norm_coord.T, axis=1)

    @property
    def volume(self):
        """Returns volume of the mesh."""
        return np.mean(self.volumes)

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        immersed_self = self.immersed_part()
        return immersed_self.surface_integral(immersed_self.faces_normals[:,2] * data, **kwargs)

    @property
    def disp_volume(self):
        return self.immersed_part().volume

    def disp_mass(self, *, rho=1000):
        return rho * self.disp_volume

    @property
    def center_of_buoyancy(self):
        """Returns center of buoyancy of the mesh."""
        immersed_self = self.immersed_part()
        coords_sq_norm = immersed_self.faces_normals * immersed_self.faces_centers**2
        return immersed_self.surface_integral(coords_sq_norm.T, axis=1) / (2*immersed_self.volume)

    @property
    def waterplane_area(self):
        """Returns water plane area of the mesh."""
        immersed_self = self.immersed_part()
        return -immersed_self.waterplane_integral(1)

    @property
    def waterplane_center(self):
        """Returns water plane center of the mesh.

        Note: Returns None if the mesh is full submerged.
        """
        immersed_self = self.immersed_part()
        waterplane_area = immersed_self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            waterplane_center = -immersed_self.waterplane_integral(
                immersed_self.faces_centers.T, axis=1) / waterplane_area
            return waterplane_center[:-1]
