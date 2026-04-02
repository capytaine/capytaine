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
        """Returns integral of given data along wet surface area.

        Parameters
        ----------
        data: np.ndarray
            Array of values at the quadrature points
            expected shape: (..., nb_faces, nb_quad_points)
        """
        return np.sum(data * self.quadrature_points[1], **kwargs)

    @property
    def wet_surface_area(self) -> float:
        """Returns wet surface area."""
        return self.immersed_part().surface_integral(1)

    @property
    def volumes(self) -> Tuple[float, float, float]:
        """Returns volumes using x, y, z components of the mesh.
        Should be the same for a regular mesh."""
        # norm_coord[i_dir, i_face, i_quad_point] = \
        #       faces_normals[i_face, i_dir] * quad_point[i_face, i_quad_point, i_dir]
        norm_coord = (
            self.faces_normals[:, None, :] * self.quadrature_points[0]
        ).transpose((2, 0, 1))
        return tuple(self.surface_integral(norm_coord, axis=(-2, -1)))

    @property
    def volume(self) -> float:
        """Returns volume of the mesh."""
        return np.mean(self.volumes)

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area.

        Parameters
        ----------
        data: np.ndarray
            Array of values at the quadrature points of the hull mesh.
            Expected shape: (..., nb_faces, nb_quad_points)
        """
        immersed_self = self.immersed_part()
        return -immersed_self.surface_integral(
            immersed_self.faces_normals[:, None, 2] * data,
            **kwargs
        )

    @property
    def disp_volume(self) -> float:
        return self.immersed_part().volume

    def disp_mass(self, *, rho=1000) -> float:
        return rho * self.disp_volume

    @property
    def center_of_buoyancy(self) -> np.ndarray:
        """Returns center of buoyancy of the mesh."""
        immersed_self = self.immersed_part()
        # coords_sq_norm[i_dir, i_face, i_quad_point] = \
        #       faces_normals[i_face, i_dir] * quad_point[i_face, i_quad_point, i_dir]**2
        coords_sq_norm = (
                immersed_self.faces_normals[:, None, :] * immersed_self.quadrature_points[0]**2
                ).transpose((2, 0, 1))
        return immersed_self.surface_integral(coords_sq_norm, axis=(-1, -2)) / (2*immersed_self.volume)

    @property
    def waterplane_area(self) -> float:
        """Returns water plane area of the mesh."""
        immersed_self = self.immersed_part()
        return immersed_self.waterplane_integral(1)

    @property
    def waterplane_center(self) -> Union[None, np.ndarray]:
        """Returns water plane center of the mesh.
        Computed as (∫x/∫1, ∫y/∫1) on the water plane.
        Returns None if the mesh is full submerged.
        """
        immersed_self = self.immersed_part()
        waterplane_area = immersed_self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            x = immersed_self.quadrature_points[0][:, :, 0]
            y = immersed_self.quadrature_points[0][:, :, 1]
            waterplane_center = (
                    immersed_self.waterplane_integral(x) / waterplane_area,
                    immersed_self.waterplane_integral(y) / waterplane_area
                    )
            return waterplane_center
