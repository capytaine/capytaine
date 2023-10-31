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

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        return self.surface_integral(self.faces_normals[:,2] * data, **kwargs)

    @property
    def wet_surface_area(self):
        """Returns wet surface area."""
        return self.surface_integral(1)

    @property
    def volumes(self):
        """Returns volumes using x, y, z components of the mesh."""
        norm_coord = self.faces_normals * self.faces_centers
        return self.surface_integral(norm_coord.T, axis=1)

    @property
    def volume(self):
        """Returns volume of the mesh."""
        return np.mean(self.volumes)

    def disp_mass(self, *, rho=1000):
        return rho * self.volume

    @property
    def center_of_buoyancy(self):
        """Returns center of buoyancy of the mesh."""
        coords_sq_norm = self.faces_normals * self.faces_centers**2
        return self.surface_integral(coords_sq_norm.T, axis=1) / (2*self.volume)

    @property
    def waterplane_area(self):
        """Returns water plane area of the mesh."""
        waterplane_area = -self.waterplane_integral(1)
        return waterplane_area

    @property
    def waterplane_center(self):
        """Returns water plane center of the mesh.

        Note: Returns None if the mesh is full submerged.
        """
        waterplane_area = self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            waterplane_center = -self.waterplane_integral(
                self.faces_centers.T, axis=1) / waterplane_area
            return waterplane_center[:-1]
