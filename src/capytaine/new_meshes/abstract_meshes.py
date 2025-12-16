from __future__ import annotations

import numpy as np

from abc import ABC, abstractmethod
from typing import Literal


class AbstractMesh(ABC):
    @property
    @abstractmethod
    def nb_vertices(self) -> int:
        ...

    @property
    @abstractmethod
    def nb_faces(self) -> int:
        ...

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

    @abstractmethod
    def merged(self) -> AbstractMesh:
        ...

    @abstractmethod
    def clipped(self, *, origin, normal, name=None) -> AbstractMesh:
        ...

    @abstractmethod
    def immersed_part(self, free_surface=0.0, *, sea_bottom=None, water_depth=None) -> AbstractMesh:
        ...
