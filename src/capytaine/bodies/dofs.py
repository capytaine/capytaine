# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


TRANSLATION_DOFS_DIRECTIONS = {"surge": (1, 0, 0), "sway": (0, 1, 0), "heave": (0, 0, 1)}
ROTATION_DOFS_AXIS = {"roll": (1, 0, 0), "pitch": (0, 1, 0), "yaw": (0, 0, 1)}


class RigidBodyDofsPlaceholder:
    """Pass an instance of this class to the FloatingBody initializer to initialize the 6 ridig body dofs."""

    def __init__(self, rotation_center=None):
        self.rotation_center = rotation_center

    def __str__(self):
        return "RigidBodyDofsPlaceholder()"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())


def rigid_body_dofs(rotation_center=None):
    return RigidBodyDofsPlaceholder(rotation_center=rotation_center)


def evaluate_translation_dof(mesh, direction=None, name=None, amplitude=1.0) -> np.array:
    if direction is None:
        if name is not None and name.lower() in TRANSLATION_DOFS_DIRECTIONS:
            direction = TRANSLATION_DOFS_DIRECTIONS[name.lower()]
        else:
            raise ValueError("A direction needs to be specified for the translation dof.")
    direction = np.asarray(direction)
    assert direction.shape == (3,)
    motion = np.empty((mesh.nb_faces, 3))
    motion[:, :] = direction
    return motion

def evaluate_rotation_dof(mesh, rotation_center=None, direction=None, name=None, amplitude=1.0) -> np.array():
    if rotation_center is None:
        rotation_center = np.array([0, 0, 0])
        LOG.warning(f"The rotation dof {name} has been initialized "
                    f"around the origin of the domain (0, 0, 0).")
    if direction is None:
        if name is not None and name.lower() in ROTATION_DOFS_AXIS:
            direction = ROTATION_DOFS_AXIS[name.lower()]
        else:
            raise ValueError("A direction needs to be specified for the rotation dof.")
    if mesh.nb_faces == 0:
        return np.empty((self.mesh.nb_faces, 3))
    else:
        motion = np.cross(rotation_center - self.mesh.faces_centers, direction)
        return amplitude * motion



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
