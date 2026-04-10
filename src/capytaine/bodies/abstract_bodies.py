"""Abstract base class for floating bodies."""
# Copyright (C) 2017-2025 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Literal

import numpy as np

from capytaine.bodies.dofs import (
    add_dofs_labels_to_vector,
    add_dofs_labels_to_matrix,
)
from capytaine.bodies.visualization import show_3d


class AbstractBody(ABC):
    """Abstract base class for FloatingBody and Multibody.

    Defines the shared interface that both single bodies and collections
    of bodies must implement.

    Subclasses must provide the following attributes (either as instance
    attributes set in __init__, or as properties/cached_properties):
        name: str
        mesh
        lid_mesh
        mesh_including_lid
        hull_mask
        dofs: dict
        mass
        center_of_mass
    """

    name: str

    # --- Dof labelling (identical in both subclasses) ---

    def add_dofs_labels_to_vector(self, vector):
        """Helper function turning a bare vector into a vector labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return add_dofs_labels_to_vector(self.dofs.keys(), vector)

    def add_dofs_labels_to_matrix(self, matrix):
        """Helper function turning a bare matrix into a matrix labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return add_dofs_labels_to_matrix(self.dofs.keys(), matrix)

    # --- Body joining (identical in both subclasses) ---

    def __add__(self, body_to_add: AbstractBody) -> AbstractBody:
        return self.join_bodies(body_to_add)

    def join_bodies(*bodies, name=None) -> AbstractBody:
        from capytaine.bodies.multibodies import Multibody
        return Multibody(bodies, name=name)

    # --- Abstract methods (different implementations) ---

    @abstractmethod
    def integrate_pressure(self, pressure): ...

    @abstractmethod
    def immersed_part(self, *args, **kwargs) -> AbstractBody: ...

    @abstractmethod
    def minimal_computable_wavelength(self): ...

    @abstractmethod
    def first_irregular_frequency_estimate(self, *args, **kwargs): ...

    @abstractmethod
    def compute_hydrostatic_stiffness(self, *, rho=1000.0, g=9.81): ...

    @abstractmethod
    def compute_rigid_body_inertia(self, rho=1000.0): ...

    @abstractmethod
    def _check_dofs_shape_consistency(self): ...

    # --- Geometric transforms ---

    @abstractmethod
    def translated(self, shift, *, name=None) -> AbstractBody: ...

    def translated_x(self, dx: float, *, name=None) -> AbstractBody:
        return self.translated([dx, 0.0, 0.0], name=name)

    def translated_y(self, dy: float, *, name=None) -> AbstractBody:
        return self.translated([0.0, dy, 0.0], name=name)

    def translated_z(self, dz: float, *, name=None) -> AbstractBody:
        return self.translated([0.0, 0.0, dz], name=name)

    @abstractmethod
    def rotated_with_matrix(self, R, *, name=None) -> AbstractBody: ...

    def rotated_x(self, angle: float, *, name=None) -> AbstractBody:
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        return self.rotated_with_matrix(R, name=name)

    def rotated_y(self, angle: float, *, name=None) -> AbstractBody:
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        return self.rotated_with_matrix(R, name=name)

    def rotated_z(self, angle: float, *, name=None) -> AbstractBody:
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        return self.rotated_with_matrix(R, name=name)

    @abstractmethod
    def mirrored(self, plane: Literal['xOz', 'yOz']) -> AbstractBody: ...

    @abstractmethod
    def clipped(self, *, origin, normal, name=None) -> AbstractBody: ...

    @abstractmethod
    def copy(self, name=None) -> AbstractBody: ...

    # --- Display ---

    @abstractmethod
    def __str__(self): ...

    @abstractmethod
    def __short_str__(self): ...

    def __repr__(self):
        return str(self)

    def show(self, *, backend=None, **kwargs):
        """Visualize the mesh using the specified backend.

        Parameters
        ----------
        backend : str, optional
            Visualization backend to use. Options are 'pyvista' or 'matplotlib'.
            By default, try several until an installed one is found.
        **kwargs
            Additional keyword arguments passed to the visualization backend.
            See :mod:`~capytaine.meshes.visualization`

        Returns
        -------
        object
            Visualization object returned by the backend (e.g., matplotlib figure).

        Raises
        ------
        NotImplementedError
            If the specified backend is not supported.
        """
        return show_3d(self, backend=backend, **kwargs)
