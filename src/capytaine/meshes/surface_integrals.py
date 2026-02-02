# Copyright 2025 Mews Labs
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

from abc import ABC
from typing import Tuple, Union

import numpy as np

class SurfaceIntegralsMixin(ABC):
    # Defines some methods inherited by AbstractMesh.
    # There are located in this other module just to make the code more tidy.

    def surface_integral(self, data, **kwargs):
        """Returns integral of given data along wet surface area."""
        return np.sum(data * self.faces_areas, **kwargs)

    @property
    def wet_surface_area(self) -> float:
        """Returns wet surface area."""
        return self.immersed_part().surface_integral(1)

    @property
    def volumes(self) -> Tuple[float, float, float]:
        """Returns volumes using x, y, z components of the mesh.
        Should be the same for a regular mesh."""
        norm_coord = self.faces_normals * self.faces_centers
        return tuple(self.surface_integral(norm_coord.T, axis=1))

    @property
    def volume(self) -> float:
        """Returns volume of the mesh."""
        return np.mean(self.volumes)

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        immersed_self = self.immersed_part()
        return immersed_self.surface_integral(immersed_self.faces_normals[:,2] * data, **kwargs)

    @property
    def disp_volume(self) -> float:
        return self.immersed_part().volume

    def disp_mass(self, *, rho=1000) -> float:
        return rho * self.disp_volume

    @property
    def center_of_buoyancy(self) -> np.ndarray:
        """Returns center of buoyancy of the mesh."""
        immersed_self = self.immersed_part()
        coords_sq_norm = immersed_self.faces_normals * immersed_self.faces_centers**2
        return immersed_self.surface_integral(coords_sq_norm.T, axis=1) / (2*immersed_self.volume)

    @property
    def waterplane_area(self) -> float:
        """Returns water plane area of the mesh."""
        immersed_self = self.immersed_part()
        return -immersed_self.waterplane_integral(1)

    @property
    def waterplane_center(self) -> Union[None, np.ndarray]:
        """Returns water plane center of the mesh.
        Returns None if the mesh is full submerged.
        """
        immersed_self = self.immersed_part()
        waterplane_area = immersed_self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            waterplane_center = -immersed_self.waterplane_integral(
                immersed_self.faces_centers.T, axis=1) / waterplane_area
            return waterplane_center[:-1]
