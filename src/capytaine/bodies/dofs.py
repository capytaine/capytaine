# Copyright 2026 Capytaine developers
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
from abc import ABC
from functools import lru_cache, partial
import logging

import numpy as np
import xarray as xr

from capytaine.meshes.abstract_meshes import AbstractMesh


LOG = logging.getLogger(__name__)


def _points_on_mesh(mesh: AbstractMesh, location: str):
    """Points of `mesh` at which the dofs are evaluated, and the index of the face containing
    each of them (an array that can be broadcast to the shape of the points
    without their last axis)."""
    # The reshapes keep the expected shapes for a mesh without faces.
    if location == "faces_centers":
        return (
            np.reshape(mesh.faces_centers, (mesh.nb_faces, 3)),
            np.arange(mesh.nb_faces)
        )
    elif location == "quadrature_points":
        return (
            np.reshape(mesh.quadrature_points[0], (*mesh.quadrature_points[1].shape, 3)),
            np.arange(mesh.nb_faces)[:, None]
        )
    elif location == "waterline":
        return (mesh.middle_edges_waterline, mesh.indices_faces_waterline)
    else:
        raise ValueError(f"Unknown location {location!r} to evaluate a dof, ",
                         "expected 'faces_centers', 'quadrature_points' or 'waterline'.")


def _overrides(dof, method_name: str) -> bool:
    """Whether the class of `dof` redefines the method `method_name` of AbstractDof."""
    return getattr(type(dof), method_name) is not getattr(AbstractDof, method_name)


class AbstractDof(ABC):
    """The representation of a degree of freedom.
    It always live within an AbstractBody, next to an AbstractMesh.

    How to use a degree of freedom
    ------------------------------
    Four magnitudes can be evaluated for a given dof:

    - the first-order motion, a vector defining how a given location of the mesh is moving,
      used for pressure integration into forces and radiation problems boundary condition;
    - its gradient (actually Jacobian matrix), used for hydrostatics, forward speed and
      second-order forces;
    - the second-order motion, the component of the body motion quadratic with respect to
      the motion amplitude, it depends on a pair of dofs and it is used for second-order
      hydrostatic forces;
    - its gradient, also used for second-order hydrostatic forces.

    The code using these classes is expected to call
    ``evaluate_[gradient_of_]motion(mesh, location)`` or
    ``evaluate_[gradient_of_]second_order_motion(other, mesh, location)``
    where ``location`` is a string which can be either ``"faces_centers"`` (the default),
    ``"quadrature_points"``, or ``"waterline"``.
    The returned array will contain the values of the fields at each chosen type of points
    of the mesh, that is of shape either ``(nb_faces, ...)``, ``(nb_faces, nb_quad_points, ...)``
    or ``(nb_edges_waterline, ...)``, where ``...`` stands for ``(3,)`` for motions
    and ``(3, 3)`` for gradients.

    How to implement a degree of freedom
    ------------------------------------
    ``CustomDof`` can be used to turn a first-order motion evaluation function in a subclass of
    ``AbstractDof``. Defining a custom dof with second-order motion requires defining a subclass
    of ``AbstractDof``, like explained below.

    For most dofs, like the built-in rigid body dofs, it is recommended to implement the low-level
    methods ``evaluate_[gradient_of_]motion_at_points`` which take an array of points of shape
    ``(..., 3)`` as inputs and return the values at these points.
    Implementing the second-order motion ``evaluate_[gradient_of_]second_order_motion_at_points``
    is optional (zero is used if none is implemented).
    These methods are meant to be private, the rest of the code is not supposed to call them,
    they might be renamed with a leading underscore in a future breaking-change version to mark that.

    Alternatively, the public API methods ``evaluate_[gradient_of_][second_order_]motion``
    can be redefined, like it is done by the ``DofOnSubmesh`` class.
    For a pair of dofs of which only one redefines ``evaluate_[gradient_of_]second_order_motion``,
    the redefined method is used whatever the order of the pair (the second-order motion is symmetric).
    For this reason, the ``*_at_points`` methods are not ``@abstractmethod`` since
    they are not strictly required.
    """

    def __str__(self):
        return repr(self)

    # ====================================== PUBLIC API =====================================
    @lru_cache
    def evaluate_motion(
        self, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        """Return the displacement as a 3 element vector for all ``location`` of ``mesh``,
        that is an array of type ``(..., 3)``."""
        points, _ = _points_on_mesh(mesh, location)
        return self.evaluate_motion_at_points(points)

    @lru_cache
    def evaluate_gradient_of_motion(
            self, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        """
        Returns a Jacobian matrix, such that output[..., i_dir, i_deriv_dir]
        is the derivative with respect to `i_deriv_dir` of the `i_dir` component of the motion.
        In other words, output[..., 0, :] is the gradient of the x-component of the motion at each
        location and output[..., :, 0] is the derivative with respect to x of the motion vector
        at each location
        """
        points, _ = _points_on_mesh(mesh, location)
        return self.evaluate_gradient_of_motion_at_points(points)

    @lru_cache
    def evaluate_second_order_motion(
            self, other: AbstractDof, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        """
        Generally a motion/deformation of a body can be expanded
        to second-order with respect to the amplitude of the dofs (q_k) as::

            δu_k = dofs[k].evaluate_motion(mesh)
            δ²u_kj = dofs[k].evaluate_second_order_motion(dofs[j], mesh)
            u = Σ_k q_k δu_k + 1/2 Σ_kj q_k q_j δ²u_kj

        ``δ²u_kj == δ²u_jk`` is the second order component coupling dof j and dof k.
        """
        if _overrides(other, "evaluate_second_order_motion"):
            # δ²u_kj is symmetric, so let `other` handle the pair with its own definition
            return other.evaluate_second_order_motion(self, mesh, location)
        points, _ = _points_on_mesh(mesh, location)
        return self.evaluate_second_order_motion_at_points(other, points)

    @lru_cache
    def evaluate_gradient_of_second_order_motion(
            self, other: AbstractDof, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        """
        Gradient (Jacobian matrix, with the same conventions as `evaluate_gradient_of_motion`)
        of δ²u_kj = dofs[k].evaluate_second_order_motion(dofs[j], mesh).
        Zero for all pair of builtin dofs, except for some pairs of rotations.
        """
        if _overrides(other, "evaluate_gradient_of_second_order_motion"):
            return other.evaluate_gradient_of_second_order_motion(self, mesh, location)
        points, _ = _points_on_mesh(mesh, location)
        return self.evaluate_gradient_of_second_order_motion_at_points(other, points)

    # =================================== LOW LEVEL HOOKS ===================================
    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        # points is an array of shape (..., 3), output should be of shape (..., 3)
        raise NotImplementedError(
            f"{type(self).__name__} does not implement `evaluate_motion_at_points`."
        )

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        # points is an array of shape (..., 3), output should be of shape (..., 3, 3)
        # See docstring of evaluate_gradient_of_motion.
        raise NotImplementedError(
            f"{type(self).__name__} does not implement `evaluate_gradient_of_motion_at_points`."
        )

    def evaluate_second_order_motion_at_points(
        self, other: AbstractDof, points: np.ndarray
    ) -> np.ndarray:
        # points is an array of shape (..., 3), output is of shape (..., 3)
        # This method returning δ²u_kj is zero by default.
        # It is overridden for some pairs of rotations.
        return np.zeros((*points.shape[:-1], 3))

    def evaluate_gradient_of_second_order_motion_at_points(
            self, other: AbstractDof, points: np.ndarray
    ) -> np.ndarray:
        # points is an array of shape (..., 3), output is of shape (..., 3, 3)
        # This method returning the gradient of δ²u_kj is zero by default.
        # It is overridden for some pairs of rotations.
        return np.zeros((*points.shape[:-1], 3, 3))


class TranslationDof(AbstractDof):
    def __init__(self, direction):
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,)

    def __repr__(self):
        return f"TranslationDof(direction={self.direction})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.tile(self.direction, (*points.shape[:-1], 1))

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.zeros((*points.shape[:-1], 3, 3))


class RotationDof(AbstractDof):
    def __init__(self, rotation_center, direction):
        self.direction = np.asarray(direction)
        assert self.direction.shape == (3,)
        if rotation_center is None:
            self.rotation_center = np.array([0.0, 0.0, 0.0])
            LOG.warning("Rigid body rotation dof has been initialized "
                        "around the origin of the domain (0, 0, 0).")
        else:
            self.rotation_center = np.asarray(rotation_center, dtype=float)
        assert self.rotation_center.shape == (3,)

    def __repr__(self):
        return f"RotationDof(rotation_center={self.rotation_center}, direction={self.direction})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        return np.cross(self.direction, points - self.rotation_center)

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        grad = np.cross(self.direction, np.eye(3)).T
        # Transposing because np.cross compute the cross product row-wise,
        # but we want it column-wise for the conventions of the Jacobian matrix.
        return np.tile(grad, (*points.shape[:-1], 1, 1))

    def _roll_pitch_yaw_rank(self) -> int:
        # 0 for roll (around x), 1 for pitch (around y), 2 for yaw (around z)
        axis = np.flatnonzero(self.direction)
        if len(axis) != 1:
            raise NotImplementedError(
                "Second order motion is only implemented for rotations around the x, y or z axes."
            )
        return int(axis[0])

    def _outer_and_inner_rotations(self, other: RotationDof):
        if not np.allclose(self.rotation_center, other.rotation_center):
            raise ValueError("For second order motion, all rotations should have the same rotation center")
        # The rotations are composed in the roll-pitch-yaw sequence of WAMIT
        # (WAMIT manual v7.6, eq. 15.48-15.49), R = R_roll R_pitch R_yaw.
        # The second-order term of the product of two different rotations depends on this order:
        # δ²u_kj = δ²u_jk = ω_outer × (ω_inner × (X - X_center))
        # with the outer rotation being the one applied last (roll, then pitch, then yaw).
        if self._roll_pitch_yaw_rank() <= other._roll_pitch_yaw_rank():
            return self, other
        else:
            return other, self

    def evaluate_second_order_motion_at_points(self, other: AbstractDof, points: np.ndarray) -> np.ndarray:
        if isinstance(other, RotationDof):
            outer, inner = self._outer_and_inner_rotations(other)
            return np.cross(outer.direction, np.cross(inner.direction, points - self.rotation_center))
        else:
            return np.zeros((*points.shape[:-1], 3))

    def evaluate_gradient_of_second_order_motion_at_points(self, other: AbstractDof, points: np.ndarray) -> np.ndarray:
        if isinstance(other, RotationDof):
            outer, inner = self._outer_and_inner_rotations(other)
            # ∇(ω_outer × (ω_inner × r)) = S(ω_outer) S(ω_inner), see `evaluate_gradient_of_motion_at_points`
            # for the transposition.
            S_outer = np.cross(outer.direction, np.eye(3)).T
            S_inner = np.cross(inner.direction, np.eye(3)).T
            return np.tile(S_outer @ S_inner, (*points.shape[:-1], 1, 1))
        else:
            return np.zeros((*points.shape[:-1], 3, 3))


class DofOnSubmesh(AbstractDof):
    """Defines a dof that is zeros everywhere except on a given set of faces
    on which another dof object is used for the definition.

    It is only defined on the mesh it has been built for, and can only be evaluated
    with the methods taking this mesh (`evaluate_motion(mesh)`, ...),
    not at arbitrary points (`evaluate_motion_at_points` is not implemented).

    Parameters
    ----------
    dof: AbstractDof
        Some other dof
    faces: boolean array
        Boolean array of size `nb_faces`,
        the total number of faces of the mesh on which this dof is defined,
        True on the faces on which `dof` should be evaluated.
    """
    def __init__(self, dof: AbstractDof, faces):
        faces = np.asarray(faces)
        if faces.dtype != bool or faces.ndim != 1:
            raise TypeError(
                "DofOnSubmesh expects a 1D boolean array of size nb_faces to define its faces, "
                f"got {faces!r}."
            )
        self.dof = dof
        self.faces = faces

    def __repr__(self):
        return f"DofOnSubmesh(dof={self.dof}, faces={self.faces})"

    def _evaluate_on_submesh(self, evaluate, value_shape, mesh, location, other=None):
        """Call `evaluate` (a method of the inner dof taking points) at the points of `location`
        on the faces of the submesh, and return 0 elsewhere. If `other` is also a `DofOnSubmesh`,
        only on the faces of both submeshes."""
        points, face_indices_of_points = _points_on_mesh(mesh, location)
        in_submesh = np.broadcast_to(self.faces[face_indices_of_points], points.shape[:-1])
        if isinstance(other, DofOnSubmesh):
            in_submesh = in_submesh & other.faces[face_indices_of_points]
        values = np.zeros((*points.shape[:-1], *value_shape))
        if np.any(in_submesh):
            # Not called on an empty submesh, where the inner dof might raise an error:
            # e.g. for two rotations of two different bodies, that have different rotation centers.
            values[in_submesh] = evaluate(points[in_submesh])
        return values

    @lru_cache
    def evaluate_motion(
            self, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        return self._evaluate_on_submesh(self.dof.evaluate_motion_at_points, (3,), mesh, location)

    @lru_cache
    def evaluate_gradient_of_motion(
            self, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        return self._evaluate_on_submesh(self.dof.evaluate_gradient_of_motion_at_points, (3, 3), mesh, location)

    @lru_cache
    def evaluate_second_order_motion(
            self, other: AbstractDof, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        inner_other = other.dof if isinstance(other, DofOnSubmesh) else other
        evaluate = partial(self.dof.evaluate_second_order_motion_at_points, inner_other)
        return self._evaluate_on_submesh(evaluate, (3,), mesh, location, other)

    @lru_cache
    def evaluate_gradient_of_second_order_motion(
            self, other: AbstractDof, mesh: AbstractMesh, location: str = "faces_centers"
    ) -> np.ndarray:
        inner_other = other.dof if isinstance(other, DofOnSubmesh) else other
        evaluate = partial(self.dof.evaluate_gradient_of_second_order_motion_at_points, inner_other)
        return self._evaluate_on_submesh(evaluate, (3, 3), mesh, location, other)


class CustomDof(AbstractDof):
    """Defines a custom dof with the same interface as the other AbstractDof.
    To be used for elastic dofs. Only first-order motion can be defined this way for now.

    Parameters
    ----------
    motion: Callable
        A function returning the motion at a given location. Required.
        The function is expected to take a 3-element vector as input (location in space)
        and return a 3-element vector (displacement at this location).
    gradient_of_motion: Callable, optional
        A function returning the gradient (actually jacobian matrix) of the motion
        The function is expected to take a 3-element vector as input (location in space)
        and return a 3×3-matrix (jacobian of the displacement at this location).
        If none is provided, zero value is used, which makes no difference for first-order hydrodynamics,
        but loses some accuracy for hydrostatics and second-order forces.
    """
    def __init__(self, motion, gradient_of_motion=None):
        self.motion = motion
        self.gradient_of_motion = gradient_of_motion

    def __repr__(self):
        grad = f", gradient_of_motion={self.gradient_of_motion}" if self.gradient_of_motion is not None else ""
        return f"CustomDof(motion={self.motion}{grad})"

    def evaluate_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        flat_points = points.reshape(-1, 3)
        motion = np.zeros((flat_points.shape[0], 3))
        for i, p in enumerate(flat_points):
            motion[i, :] = self.motion(p)
        return motion.reshape(*points.shape[:-1], 3)

    def evaluate_gradient_of_motion_at_points(self, points: np.ndarray) -> np.ndarray:
        flat_points = points.reshape(-1, 3)
        grad = np.zeros((flat_points.shape[0], 3, 3))
        if self.gradient_of_motion is not None:
            for i, p in enumerate(flat_points):
                grad[i, :, :] = self.gradient_of_motion(p)
        return grad.reshape(*points.shape[:-1], 3, 3)


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
