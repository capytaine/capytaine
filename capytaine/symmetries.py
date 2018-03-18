#!/usr/bin/env python
# coding: utf-8
"""Special bodies using their symmetries to speed up the computations.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.geometry import Plane

from capytaine.bodies import FloatingBody
from capytaine.bodies_collection import CollectionOfFloatingBodies


LOG = logging.getLogger(__name__)


class _SymmetricBody(CollectionOfFloatingBodies):
    def __add__(self, body_to_add):
        return CollectionOfFloatingBodies([self, body_to_add])


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)


class ReflectionSymmetry(_SymmetricBody):
    """A body composed of two symmetrical halves."""

    def __init__(self, half, plane, name=None):
        """Initialize the body.

        Parameters
        ----------
        half : FloatingBody
            a FloatingBody instance describing half of the body
        plane : Plane
            the symmetry plane across which the half body is mirrored
        name : string, optional
            a name for the body
        """
        assert isinstance(half, FloatingBody)
        assert isinstance(plane, Plane)
        assert plane.normal[2] == 0  # Only vertical reflection planes are supported

        half.nb_matrices_to_keep *= 2

        self.plane = plane

        other_half = half.copy()
        other_half.mirror(plane)
        other_half.name = "mirror_of_" + half.name

        CollectionOfFloatingBodies.__init__(self, [half, other_half])

        if name is None:
            self.name = f"ReflectionSymmetry({half.name})"
        else:
            self.name = name
        LOG.info(f"New mirror symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in half.dofs.items():
            self.dofs['mirrored_' + name] = np.concatenate([dof, dof])

    def get_immersed_part(self, **kwargs):
        return ReflectionSymmetry(self.subbodies[0].get_immersed_part(**kwargs),
                                  plane=self.plane,
                                  name=f"{self.name}_clipped")


class TranslationalSymmetry(_SymmetricBody):
    """A body composed of a pattern repeated and translated."""

    def __init__(self, body_slice, translation, nb_repetitions=1, name=None):
        """Initialize the body.

        Parameters
        ----------
        body_slice : FloatingBody
            the pattern that will be repeated to form the whole body
        translation : array(3)
            the vector of the translation
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : string, optional
            a name for the body
        """
        assert isinstance(body_slice, FloatingBody)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        translation = np.asarray(translation)
        assert translation.shape == (3,)
        assert translation[2] == 0  # Only horizontal translation are supported.

        self.translation = translation

        body_slice.nb_matrices_to_keep *= nb_repetitions+1
        slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy(name=f"repetition_{i}_of_{body_slice.name}")
            new_slice.translate(i*translation)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            slices.append(new_slice)

        CollectionOfFloatingBodies.__init__(self, slices)

        if name is None:
            self.name = f"TranslationSymmetry({body_slice.name})"
        else:
            self.name = name
        LOG.info(f"New translation symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["translated_" + name] = np.concatenate([dof]*nb_repetitions)

    def get_immersed_part(self, **kwargs):
        return TranslationalSymmetry(self.subbodies[0].get_immersed_part(**kwargs),
                                     translation=self.translation,
                                     nb_repetitions=self.nb_subbodies-1,
                                     name=f"{self.name}_clipped")


class AxialSymmetry(_SymmetricBody):
    """A body composed of a pattern rotated around a vertical axis."""

    def __init__(self, body_slice, point_on_rotation_axis=np.zeros(3), nb_repetitions=1, name=None):
        """Initialize the body.

        Parameters
        ----------
        body_slice : FloatingBody
            the pattern that will be repeated to form the whole body
        point_on_rotation_axis : array(3)
            one point on the rotation axis. The axis is supposed to be vertical.
            TODO: Use an Axis class.
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : string, optional
            a name for the body
        """
        assert isinstance(body_slice, FloatingBody)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        point_on_rotation_axis = np.asarray(point_on_rotation_axis)
        assert point_on_rotation_axis.shape == (3,)
        self.point_on_rotation_axis = point_on_rotation_axis

        body_slice.nb_matrices_to_keep *= nb_repetitions+1
        slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy(name=f"rotation_{i}_of_{body_slice.name}")
            new_slice.translate(-point_on_rotation_axis)
            new_slice.rotate_z(2*i*np.pi/(nb_repetitions+1))
            new_slice.translate(point_on_rotation_axis)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            slices.append(new_slice)

        CollectionOfFloatingBodies.__init__(self, slices)

        if name is None:
            self.name = f"AxialSymmetry({body_slice.name})"
        else:
            self.name = name
        LOG.info(f"New rotation symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["rotated_" + name] = np.concatenate([dof]*nb_repetitions)

    @staticmethod
    def from_profile(profile,
                     z_range=np.linspace(-5, 0, 20),
                     point_on_rotation_axis=np.zeros(3),
                     nphi=20,
                     name=None):
        """Return a floating body using the axial symmetry.
        The shape of the body can be defined either with a function defining the profile as [f(z), 0, z] for z in z_range.
        Alternatively, the profile can be defined as a list of points.
        The number of vertices along the vertical direction is len(z_range) in the first case and profile.shape[0] in the second case.

        Parameters
        ----------
        profile : function(float) → float  or  array(N, 3)
            define the shape of the body either as a function or a list of points.
        z_range: array(N)
            used only if the profile is defined as a function.
        point_on_rotation_axis: array(3)
            a single point to define the rotation axis (the direction is always vertical)
        nphi : int
            number of vertical slices forming the body
        name : string
            a name identifying the body (default: "repeated_slice_id" where id is an unique integer).

        Returns
        -------
        AxialSymmetry
            the generated body
        """
        if name is None:
            name = f"axi-symmetric_body_{next(Mesh._ids)}"

        if callable(profile):
            x_values = [profile(z) for z in z_range]
            profile_array = np.stack([x_values, np.zeros(len(z_range)), z_range]).T
        else:
            profile_array = np.asarray(profile)
        assert len(profile_array.shape) == 2
        assert profile_array.shape[1] == 3

        n = profile_array.shape[0]
        angle = 2 * np.pi / nphi

        rotated_profile = FloatingBody(Mesh(profile_array, np.zeros((0, 4)), name="rotated_profile_mesh"), name="rotated_profile")
        rotated_profile.rotate_z(angle)

        nodes_slice = np.concatenate([profile_array, rotated_profile.vertices])
        faces_slice = np.array([[i, i+n, i+n+1, i+1] for i in range(n-1)])
        body_slice = FloatingBody(Mesh(nodes_slice, faces_slice, name=f"slice_of_{name}_mesh"), name=f"slice_of_{name}")
        body_slice.mesh.merge_duplicates()
        body_slice.mesh.heal_triangles()

        return AxialSymmetry(body_slice, point_on_rotation_axis=point_on_rotation_axis, nb_repetitions=nphi-1, name=name)

    def get_immersed_part(self, **kwargs):
        return AxialSymmetry(self.subbodies[0].get_immersed_part(**kwargs),
                             point_on_rotation_axis=self.point_on_rotation_axis,
                             nb_repetitions=self.nb_subbodies-1,
                             name=f"{self.name}_clipped")

