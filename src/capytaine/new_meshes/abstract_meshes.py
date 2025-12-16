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

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Literal

import numpy as np

from capytaine.tools.deprecation_handling import _get_water_depth

LOG = logging.getLogger(__name__)

class AbstractMesh(ABC):
    @property
    @abstractmethod
    def nb_vertices(self) -> int:
        ...

    @property
    @abstractmethod
    def nb_faces(self) -> int:
        ...

    @property
    @abstractmethod
    def faces_normals(self) -> np.ndarray:
        ...

    @property
    @abstractmethod
    def faces_areas(self) -> np.ndarray:
        ...

    @property
    @abstractmethod
    def faces_centers(self) -> np.ndarray:
        ...

    @property
    @abstractmethod
    def faces_radiuses(self) -> np.ndarray:
        ...

    @property
    @abstractmethod
    def faces(self) -> np.ndarray:
        ...

    @property
    @abstractmethod
    def quadrature_points(self) -> np.ndarray:
        ...

    @cached_property
    def z_span(self) -> Tuple[float, float]:
        return (self.vertices[:, 2].min(), self.vertices[:, 2].max())

    @abstractmethod
    def __str__(self) -> str:
        ...

    @abstractmethod
    def __short_str__(self) -> str:
        ...

    @abstractmethod
    def extract_faces(self, faces_id, *, name=None) -> AbstractMesh:
        ...

    @abstractmethod
    def translated(self, shift, *, name=None) -> AbstractMesh:
        ...

    def translated_x(self, dx: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh translated in the x-direction along `dx`."""
        return self.translated([dx, 0.0, 0.0], name=name)

    def translated_y(self, dy: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh translated in the y-direction along `dy`."""
        return self.translated([0.0, dy, 0.0], name=name)

    def translated_z(self, dz: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh translated in the z-direction along `dz`."""
        return self.translated([0.0, 0.0, dz], name=name)

    @abstractmethod
    def rotated_with_matrix(self, R, *, name=None) -> AbstractMesh:
        ...

    def rotated_x(self, angle: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh rotated around the x-axis using the provided rotation angle in radians"""
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        return self.rotated_with_matrix(R, name=name)

    def rotated_y(self, angle: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh rotated around the y-axis using the provided rotation angle in radians"""
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        return self.rotated_with_matrix(R, name=name)

    def rotated_z(self, angle: float, *, name=None) -> AbstractMesh:
        """Return a new Mesh rotated around the z-axis using the provided rotation angle in radians"""
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        return self.rotated_with_matrix(R, name=name)

    def mirrored(self, plane: Literal['xOz', 'yOz'], *, name=None) -> AbstracMesh:
        ...

    @abstractmethod
    def join_meshes(*meshes, return_masks=False, name=None) -> AbstractMesh:
        ...

    def _common_metadata_keys(*meshes):
        metadata_keys = [set(m.faces_metadata.keys()) for m in meshes]
        common_metadata_keys = set.intersection(*metadata_keys)
        lost_metadata_keys = set.union(*metadata_keys) - common_metadata_keys
        if len(lost_metadata_keys) > 0:
            LOG.warning(f'The following metadata have been dropped when joining meshes: {lost_metadata_keys}')
        return common_metadata_keys

    def __add__(self, other: AbstractMesh) -> AbstractMesh:
        """Combine two meshes using the + operator.

        Parameters
        ----------
        other : Mesh
            Another mesh to combine with this one.

        Returns
        -------
        Mesh
            New mesh containing vertices and faces from both meshes.
        """
        if self.name is not None or other.name is not None:
            name = f"{self.name}+{other.name}"
        else:
            name = None
        return self.join_meshes(other, name=name)

    @abstractmethod
    def with_normal_vector_going_down(self, **kwargs) -> AbstractMesh:
        ...

    @abstractmethod
    def copy(self) -> AbstractMesh:
        ...

    def with_metadata(self, **new_metadata) -> AbstractMesh:
        faces_metadata = self.faces_metadata.copy()
        for k, v in new_metadata.items():
            faces_metadata[k] = v
        return self.copy(faces_metadata=faces_metadata)

    def pop_metadata(self, metadata_name) -> Tuple[AbstractMesh, np.ndarray]:
        faces_metadata = self.faces_metadata.copy()
        data = faces_metadata.pop(metadata_name)
        return self.copy(faces_metadata=faces_metadata), data

    def without_metadata(self, *metadata_names) -> AbstractMesh:
        faces_metadata = self.faces_metadata.copy()
        for k in metadata_names:
            del faces_metadata[k]
        return self.copy(faces_metadata=faces_metadata)

    def without_any_metadata(self) -> AbstractMesh:
        return self.copy(faces_metadata={})

    @abstractmethod
    def merged(self) -> AbstractMesh:
        ...

    @abstractmethod
    def clipped(self, *, origin, normal, name=None) -> AbstractMesh:
        ...

    def immersed_part(self, free_surface=0.0, *, sea_bottom=None, water_depth=None) -> AbstractMesh:
        """
        Clip the mesh to keep only the part below the free surface.

        Parameters
        ----------
        free_surface: float
            The :math:`z` coordinate of the free surface (default: 0.0)
        water_depth: Optional[float]
            The water depth, as a positive value (default: infinity)

        Returns
        -------
        Mesh
            A new Mesh instance that has been clipped.
        """
        water_depth = _get_water_depth(free_surface, water_depth, sea_bottom,
                                       default_water_depth=np.inf)
        if (free_surface - water_depth <= self.z_span[0]
            and self.z_span[1] <= free_surface):  # Already clipped
            return self  # Shortcut for performance
        clipped = self.clipped(origin=(0, 0, 0), normal=(0, 0, 1))
        if water_depth < np.inf:
            clipped = clipped.clipped(origin=(0, 0, free_surface-water_depth), normal=(0, 0, -1))
        return clipped

    @abstractmethod
    def show(self, *, backend=None, **kwargs):
        ...

    def show_pyvista(self, **kwargs):
        """
        Equivalent to show(backend="pyvista").
        See also :func:`~capytaine.new_meshes.visualization.show_pyvista`
        """
        return self.show(backend="pyvista", **kwargs)

    def show_matplotlib(self, **kwargs):
        """
        Equivalent to show(backend="matplotlib").
        See also :func:`~capytaine.new_meshes.visualization.show_matplotlib`
        """
        return self.show(backend="matplotlib", **kwargs)
