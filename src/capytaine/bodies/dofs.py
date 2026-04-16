# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC, abstractmethod
from functools import lru_cache
import logging

import numpy as np
import xarray as xr

from capytaine.meshes.abstract_meshes import AbstractMesh

LOG = logging.getLogger(__name__)


class AbstractDof(ABC):
    @lru_cache
    def evaluate_motion(self, mesh: AbstractMesh) -> np.ndarray:
        if mesh.nb_faces == 0:
            return np.empty((mesh.nb_faces, 3))
        else:
            return self.evaluate_motion_at_points(mesh.faces_centers)

    @abstractmethod
    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        # points is an array of shape (nb_points, 3)
        # output is of shape (nb_points, 3)
        ...

    @lru_cache
    def evaluate_gradient_of_motion(self, mesh: AbstractMesh) -> np.ndarray:
        # output is of shape (nb_points, 3, 3)
        if mesh.nb_faces == 0:
            return np.empty((mesh.nb_faces, 3, 3))
        else:
            return self.evaluate_gradient_of_motion_at_points(mesh.faces_centers)

    @abstractmethod
    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        # points is an array of shape (nb_points, 3)
        # output is of shape (nb_points, 3, 3)
        # output is a Jacobian matrix, such that output[i_point, i_dir, i_deriv_dir]
        # is the derivative with respect to `i_deriv_dir` of the `i_dir` component of the motion.
        # In other words, output[:, 0, :] is the gradient of the x-component of the motion on each face and
        # output[:, :, 0] is the derivative with respect to x of the motion vector on each face.
        ...


class TranslationDof(AbstractDof):
    def __init__(self, direction):
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,)

    def __str__(self):
        return f"TranslationDof(direction={self.direction})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.tile(self.direction, (points.shape[0], 1))

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.zeros((points.shape[0], 3, 3))


class RotationDof(AbstractDof):
    def __init__(self, rotation_center, direction):
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,)
        if rotation_center is None:
            self.rotation_center = np.array([0, 0, 0])
            LOG.warning("Rigid body rotation dof has been initialized "
                        "around the origin of the domain (0, 0, 0).")
        else:
            self.rotation_center = np.asarray(rotation_center)
        assert self.rotation_center.shape == (3,)

    def __str__(self):
        return f"RotationDof(rotation_center={self.rotation_center}, direction={self.direction})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.cross(self.direction, points - self.rotation_center)

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        grad = np.cross(self.direction, np.eye(3)).T
        # Transposing because np.cross compute the cross product row-wise,
        # but we want it column-wise for the conventions of the Jacobian matrix.
        return np.tile(grad, (points.shape[0], 1, 1))


class DofOnSubmesh(AbstractDof):
    """Defines a dof that is zeros everywhere except on a given set of faces
    on which another dof object is used for the definition.

    Parameters
    ----------
    dof: AbstractDof
        Some other dof
    faces: boolean array (or slice or range)
        The indices of the faces on which the `dof` should be evaluated.
        Can be provided as a list of indices or a slice,
        or as a boolean array of size `nb_faces` the total size of the mesh on which this dof is defined.
    """
    def __init__(self, dof: AbstractDof, faces):
        self.dof = dof
        self.faces = faces

    def __str__(self):
        return f"DofOnSubmesh(dof={self.dof}, faces={self.faces})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        motion = np.zeros((points.shape[0], 3))
        motion[self.faces, :] = self.dof.evaluate_motion_at_points(points[self.faces, :])
        return motion

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        grad = np.zeros((points.shape[0], 3, 3))
        grad[self.faces, :, :] = self.dof.evaluate_gradient_of_motion_at_points(points[self.faces, :])
        return grad


def is_rigid_body_dof(dof):
    return (
            isinstance(dof, TranslationDof)
            or isinstance(dof, RotationDof)
            # or (isinstance(dof, DofOnSubmesh) and is_rigid_body_dof(dof.dof))
            )


def rigid_body_dofs(only=None, rotation_center=None):
    """Pass this to FloatingBody initializer to give it rigid body dofs.

    Parameters
    ----------
    only: sequence of str, optional
        list of the name of the rigid body dofs to be included.
        By default: all six of them
    rotation_center: np.array, optional
        the center for the definition of the rotations
    """
    if rotation_center is None:
        rotation_center = np.array([0, 0, 0])

        if not (only is not None
                and set(only).issubset({"Surge", "Sway", "Heave"})
                ): # Skip the warning if only translations are required.
            LOG.warning("Rigid body rotation dofs have been initialized "
                        "around the origin of the domain (0, 0, 0).")
        # This warning is redundant with the one in RotationDof.__init__,
        # but it is done here to have a single warning displayed on screen
        # when a rigid body is initialized.

    dofs = {
        "Surge": TranslationDof(direction=(1, 0, 0)),
        "Sway": TranslationDof(direction=(0, 1, 0)),
        "Heave": TranslationDof(direction=(0, 0, 1)),
        "Roll": RotationDof(rotation_center=rotation_center, direction=(1, 0, 0)),
        "Pitch": RotationDof(rotation_center=rotation_center, direction=(0, 1, 0)),
        "Yaw": RotationDof(rotation_center=rotation_center, direction=(0, 0, 1)),
        }
    if only is not None:
        dofs = {k: v for k, v in dofs.items() if k in only}
    return dofs


def normalize_name(name):
    return name[0].upper() + name[1:].lower()


def add_dofs_labels_to_vector(dof_names, vector):
    """Helper function turning a bare vector into a vector labelled by the name of the dofs,
    to be used for instance for the computation of RAO."""
    return xr.DataArray(data=np.asarray(vector), dims=['influenced_dof'],
                        coords={'influenced_dof': list(dof_names)},
                        )

def add_dofs_labels_to_matrix(dof_names, matrix):
    """Helper function turning a bare matrix into a matrix labelled by the name of the dofs,
    to be used for instance for the computation of RAO."""
    return xr.DataArray(data=np.asarray(matrix), dims=['influenced_dof', 'radiating_dof'],
                        coords={'influenced_dof': list(dof_names), 'radiating_dof': list(dof_names)},
                        )
