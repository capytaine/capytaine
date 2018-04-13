#!/usr/bin/env python
# coding: utf-8
"""Special meshes with symmetries, useful to speed up the computations.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.geometry import Plane

from capytaine.meshes_collection import CollectionOfMeshes

LOG = logging.getLogger(__name__)


class SymmetricMesh(CollectionOfMeshes):
    pass


class ReflectionSymmetry(SymmetricMesh):
    def __init__(self, half, plane, name=None):
        """A mesh with one vertical symmetry plane.

        Parameters
        ----------
        half : Mesh or CollectionOfMeshes
            a mesh describing half of the body
        plane : Plane
            the symmetry plane across which the half body is mirrored
        name :str, optional
            a name for the mesh
        """
        assert isinstance(half, Mesh) or isinstance(half, CollectionOfMeshes)
        assert isinstance(plane, Plane)
        assert plane.normal[2] == 0  # Only vertical reflection planes are supported

        self.plane = plane

        other_half = half.copy()
        other_half.mirror(plane)
        other_half.name = "mirror_of_" + half.name

        CollectionOfMeshes.__init__(self, (half, other_half))

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, half.name)
        else:
            self.name = name
        LOG.info(f"New mirror symmetric mesh: {self.name}.")

    # def get_clipped_mesh(self, **kwargs):
    #     return ReflectionSymmetry(self.subbodies[0].get_clipped_mesh(**kwargs),
    #                               plane=self.plane,
    #                               name=f"{self.name}_clipped")



class TranslationalSymmetry(SymmetricMesh):
    def __init__(self, mesh_slice, translation, nb_repetitions=1, name=None):
        """A mesh with a repeating pattern by translation.

        Parameters
        ----------
        mesh_slice : Mesh or CollectionOfMeshes
            the pattern that will be repeated to form the whole body
        translation : array(3)
            the vector of the translation
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : str, optional
            a name for the mesh
        """
        assert isinstance(mesh_slice, Mesh) or isinstance(mesh_slice, CollectionOfMeshes)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        translation = np.asarray(translation)
        assert translation.shape == (3,)
        assert translation[2] == 0  # Only horizontal translation are supported.

        self.translation = translation

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = mesh_slice.copy()
            new_slice.name = f"repetition_{i}_of_{mesh_slice.name}"
            new_slice.translate(i*translation)
            slices.append(new_slice)

        CollectionOfMeshes.__init__(self, slices)

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, mesh_slice.name)
        else:
            self.name = name
        LOG.info(f"New translation symmetric mesh: {self.name}.")

    # def get_clipped_mesh(self, **kwargs):
    #     return TranslationalSymmetry(self.subbodies[0].get_clipped_mesh(**kwargs),
    #                                  translation=self.translation,
    #                                  nb_repetitions=self.nb_subbodies-1,
    #                                  name=f"{self.name}_clipped")


class AxialSymmetry(SymmetricMesh):
    def __init__(self, mesh_slice, point_on_rotation_axis=np.zeros(3), nb_repetitions=1, name=None):
        """A mesh with a repeating pattern by rotation.

        Parameters
        ----------
        mesh_slice : Mesh or CollectionOfMeshes
            the pattern that will be repeated to form the whole body
        point_on_rotation_axis : array(3)
            one point on the rotation axis. The axis is supposed to be vertical.
            TODO: Use an Axis class.
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : str, optional
            a name for the mesh
        """
        assert isinstance(mesh_slice, Mesh) or isinstance(mesh_slice, CollectionOfMeshes)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        point_on_rotation_axis = np.asarray(point_on_rotation_axis)
        assert point_on_rotation_axis.shape == (3,)

        self.point_on_rotation_axis = point_on_rotation_axis

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = mesh_slice.copy()
            new_slice.name = f"rotation_{i}_of_{mesh_slice.name}"
            new_slice.translate(-point_on_rotation_axis)
            new_slice.rotate_z(2*i*np.pi/(nb_repetitions+1))
            new_slice.translate(point_on_rotation_axis)
            slices.append(new_slice)

        CollectionOfMeshes.__init__(self, tuple(slices))

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, mesh_slice.name)
        else:
            self.name = name
        LOG.info(f"New rotation symmetric mesh: {self.name}.")

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
        profile : function(float → float)  or  array(N, 3)
            define the shape of the body either as a function or a list of points.
        z_range: array(N), optional
            used only if the profile is defined as a function.
        point_on_rotation_axis: array(3), optional
            a single point to define the rotation axis (the direction is always vertical)
        nphi : int, optional
            number of vertical slices forming the body
        name : str, optional
            name of the generated body (optional)

        Returns
        -------
        AxialSymmetry
            the generated mesh
        """

        if name is None:
            name = "axisymmetric_mesh"

        if callable(profile):
            x_values = [profile(z) for z in z_range]
            profile_array = np.stack([x_values, np.zeros(len(z_range)), z_range]).T
        else:
            profile_array = np.asarray(profile)

        assert len(profile_array.shape) == 2
        assert profile_array.shape[1] == 3

        n = profile_array.shape[0]
        angle = 2 * np.pi / nphi

        rotated_profile = Mesh(profile_array, np.zeros((0, 4)), name="rotated_profile_mesh")
        rotated_profile.rotate_z(angle)

        nodes_slice = np.concatenate([profile_array, rotated_profile.vertices])
        faces_slice = np.array([[i, i+n, i+n+1, i+1] for i in range(n-1)])
        body_slice = Mesh(nodes_slice, faces_slice, name=f"slice_of_{name}_mesh")
        body_slice.merge_duplicates()
        body_slice.heal_triangles()

        return AxialSymmetry(body_slice, point_on_rotation_axis=point_on_rotation_axis, nb_repetitions=nphi-1, name=name)

    # def get_clipped_mesh(self, **kwargs):
    #     return AxialSymmetry(self.subbodies[0].get_clipped_mesh(**kwargs),
    #                          point_on_rotation_axis=self.point_on_rotation_axis,
    #                          nb_repetitions=self.nb_subbodies-1,
    #                          name=f"{self.name}_clipped")
