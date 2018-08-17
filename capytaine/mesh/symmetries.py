#!/usr/bin/env python
# coding: utf-8
"""Special meshes with symmetries, useful to speed up the computations."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging
from functools import wraps

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.mesh.meshes_collection import CollectionOfMeshes
from capytaine.tools.geometry import Plane

LOG = logging.getLogger(__name__)


class SymmetricMesh(CollectionOfMeshes):
    pass


class ReflectionSymmetry(SymmetricMesh):
    def __new__(cls, half, plane, name=None):
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

        other_half = half.copy()
        other_half.mirror(plane)
        other_half.name = "mirror_of_" + half.name

        self = super().__new__(cls, (half, other_half))

        self.plane = plane

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, half.name)
        else:
            self.name = name

        LOG.info(f"New mirror symmetric mesh: {self.name}.")

        return self

    @property
    def half(self):
        return self[0]

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.name + '\n' + ' ├─' + self.half.tree_view().replace('\n', '\n | ') + '\n'
                    + f" └─mirrored copy of the above {self.half.name}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return ReflectionSymmetry(self.half.copy(), self.plane, name=self.name)

    def get_immersed_part(self, **kwargs):
        clipped_submesh = self.half.get_immersed_part(self, **kwargs)
        if clipped_submesh is not None:
            return ReflectionSymmetry(clipped_submesh,
                                      plane=self.plane,
                                      name=f"{self.name}_clipped")
        else:
            return None


class TranslationalSymmetry(SymmetricMesh):
    def __new__(cls, mesh_slice, translation, nb_repetitions=1, name=None):
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

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = mesh_slice.copy()
            new_slice.name = f"repetition_{i}_of_{mesh_slice.name}"
            new_slice.translate(i*translation)
            slices.append(new_slice)

        self = super().__new__(cls, slices)

        self.translation = translation

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, mesh_slice.name)
        else:
            self.name = name
        LOG.info(f"New translation symmetric mesh: {self.name}.")

        return self

    @property
    def first_slice(self):
        return self[0]

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.name + '\n' + ' ├─' + self.first_slice.tree_view().replace('\n', '\n | ') + '\n'
                    + f" └─{len(self)-1} translated copies of the above {self.first_slice.name}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return TranslationalSymmetry(self.first_slice.copy(), self.translation, nb_repetitions=len(self)-1, name=self.name)

    def get_immersed_part(self, **kwargs):
        clipped_submesh = self.first_slice.get_immersed_part(**kwargs)
        if clipped_submesh is not None:
            return TranslationalSymmetry(clipped_submesh,
                                         translation=self.translation,
                                         nb_repetitions=self.nb_submeshes-1,
                                         name=f"{self.name}_clipped")
        else:
            return None

    def join(*list_of_symmetric_meshes):
        """Experimental routine to merge similar symmetries."""
        assert all([isinstance(mesh, TranslationalSymmetry) for mesh in list_of_symmetric_meshes])
        assert all([np.allclose(list_of_symmetric_meshes[0].translation, mesh.translation) for mesh in list_of_symmetric_meshes[1:]])
        assert all([list_of_symmetric_meshes[0].nb_submeshes == mesh.nb_submeshes for mesh in list_of_symmetric_meshes[1:]])
        return TranslationalSymmetry(
            sum([mesh.first_slice.merge() for mesh in list_of_symmetric_meshes[1:]],
                list_of_symmetric_meshes[0].first_slice.merge()),
            list_of_symmetric_meshes[0].translation, list_of_symmetric_meshes[0].nb_submeshes,
        )


class AxialSymmetry(SymmetricMesh):
    def __new__(cls, mesh_slice, point_on_rotation_axis=np.zeros(3), nb_repetitions=1, name=None):
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

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = mesh_slice.copy()
            new_slice.name = f"rotation_{i}_of_{mesh_slice.name}"
            new_slice.translate(-point_on_rotation_axis)
            new_slice.rotate_z(2*i*np.pi/(nb_repetitions+1))
            new_slice.translate(point_on_rotation_axis)
            slices.append(new_slice)

        self = super().__new__(cls, slices)

        self.point_on_rotation_axis = point_on_rotation_axis

        if name is None:
            self.name = CollectionOfMeshes.format_name(self, mesh_slice.name)
        else:
            self.name = name
        LOG.info(f"New rotation symmetric mesh: {self.name}.")

        return self

    @property
    def first_slice(self):
        return self[0]

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.name + '\n' + ' ├─' + self.first_slice.tree_view().replace('\n', '\n | ') + '\n'
                    + f" └─{len(self)-1} rotated copies of the above {self.first_slice.name}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return AxialSymmetry(self.first_slice.copy(), self.point_on_rotation_axis, nb_repetitions=len(self)-1, name=self.name)

    def get_immersed_part(self, **kwargs):
        clipped_submesh = self.first_slice.get_immersed_part(**kwargs)
        if clipped_submesh is not None:
            return AxialSymmetry(clipped_submesh,
                                 point_on_rotation_axis=self.point_on_rotation_axis,
                                 nb_repetitions=self.nb_submeshes-1,
                                 name=f"{self.name}_clipped")
        else:
            return None

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


###################################################################################################
#                          Application for hierarchical matrix building                           #
###################################################################################################

def use_symmetries(build_matrices):
    """Decorator for the matrix building functions.

    Parameters
    ----------
    build_matrices: function
        Function that takes as argument two meshes and several other parameters and returns an
        influence matrix.

    Returns
    -------
    function
        A similar function that returns a block matrix based on the symmetries of the meshes.
    """
    from capytaine.Toeplitz_matrices import BlockCirculantMatrix, BlockToeplitzMatrix

    @wraps(build_matrices)
    def build_matrices_with_symmetries(solver, mesh1, mesh2, *args, _rec_depth=1):
        """Assemble the influence matrices using symmetries of the body.⎈

        The method is basically an ugly multiple dispatch on the kind of bodies.
        For symmetric structures, the method is called recursively on all of the sub-bodies.

        Parameters
        ----------
        solver
            Passed to the actual evaluation of the coefficients
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        *args
            Passed to the actual evaluation of the coefficients
        _rec_depth: tuple, optional
            internal parameter: recursion accumulator for pretty log printing

        Returns
        -------
        matrix-like
            influence matrix (integral of the Green function)
        """

        if (isinstance(mesh1, ReflectionSymmetry)
                and isinstance(mesh2, ReflectionSymmetry)
                and mesh1.plane == mesh2.plane):

            LOG.debug("\t" * _rec_depth +
                      f"Evaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      f"using mirror symmetry.")

            S_a, V_a = build_matrices_with_symmetries(
                solver, mesh1[0], mesh2[0], *args,
                _rec_depth=_rec_depth+1)
            S_b, V_b = build_matrices_with_symmetries(
                solver, mesh1[0], mesh2[1], *args,
                _rec_depth=_rec_depth+1)

            return BlockToeplitzMatrix([S_a, S_b]), BlockToeplitzMatrix([V_a, V_b])

        elif (isinstance(mesh1, TranslationalSymmetry)
              and isinstance(mesh2, TranslationalSymmetry)
              and np.allclose(mesh1.translation, mesh2.translation)
              and mesh1.nb_submeshes == mesh2.nb_submeshes):

            LOG.debug("\t" * _rec_depth +
                      f"Evaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                      "using translational symmetry.")

            S_list, V_list = [], []
            for submesh in mesh2:
                S, V = build_matrices_with_symmetries(
                    solver, mesh1[0], submesh, *args,
                    _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)

            return BlockToeplitzMatrix(S_list), BlockToeplitzMatrix(V_list)

        elif (isinstance(mesh1, AxialSymmetry)
              and mesh1 is mesh2):  # TODO: Generalize: if mesh1.axis == mesh2.axis

            LOG.debug("\t" * _rec_depth +
                      f"Evaluating matrix of {mesh1.name} on itself "
                      f"using rotation symmetry.")

            S_list, V_list = [], []
            for submesh in mesh2[:mesh2.nb_submeshes // 2 + 1]:
                S, V = build_matrices_with_symmetries(
                    solver, mesh1[0], submesh, *args,
                    _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)

            return BlockCirculantMatrix(S_list, size=mesh1.nb_submeshes), BlockCirculantMatrix(V_list, size=mesh1.nb_submeshes)

        else:
            #Actual evaluation of coefficients using the Green function.
            LOG.debug("\t" * _rec_depth +
                      f"Computation of the Green function between {mesh1.name} "
                      f"and {'itself' if mesh2 is mesh1 else mesh2.name}.")

            return build_matrices(solver, mesh1, mesh2, *args)

    return build_matrices_with_symmetries

