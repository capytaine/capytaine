# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

from abc import ABC
import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


class AbstractDof(ABC):
    ...


class TranslationDof(AbstractDof):
    def __init__(self, direction, amplitude=1.0):
        self.direction = np.asarray(direction)
        assert direction.shape == (3,)
        self.amplitude = amplitude

    def evaluate_motion(self, mesh) -> np.array:
        motion = np.empty((mesh.nb_faces, 3))
        motion[:, :] = self.direction
        return self.amplitude * motion


class RotationDof(AbstractDof):
    def __init__(self, rotation_center, direction, amplitude=1.0):
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,)
        self.rotation_center = np.asarray(rotation_center)
        assert self.rotation_center.shape == (3,)
        self.amplitude = amplitude

    def evaluate_motion(self, mesh) -> np.array:
        if mesh.nb_faces == 0:
            return np.empty((mesh.nb_faces, 3))
        else:
            motion = np.cross(self.rotation_center - mesh.faces_centers, self.direction)
            return self.amplitude * motion


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
        LOG.warning("Rigid body rotation dofs have been initialized "
                    "around the origin of the domain (0, 0, 0).")
    dofs = {
        "Surge": TranslationDof(direction=(1, 0, 0)),
        "Sway": TranslationDof(direction=(0, 1, 0)),
        "Heave": TranslationDof(direction=(0, 0, 1)),
        "Roll": RotationDof(rotation_center=rotation_center, direction=(1, 0, 0)),
        "Pitch": RotationDof(rotation_center=rotation_center, direction=(0, 1, 0)),
        "Yaw": RotationDof(rotation_center=rotation_center, direction=(0, 1, 0)),
        }
    if only is not None:
        dofs = {k: v for k, v in dofs.items() if k in only}
    return dofs


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
