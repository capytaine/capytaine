#!/usr/bin/env python
# coding: utf-8
"""Storing a set of meshes."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

from itertools import chain, accumulate

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.tools.geometry import Abstract3DObject, inplace_transformation

NAME_MAX_LENGTH = 180


class CollectionOfMeshes(tuple, Abstract3DObject):
    """A tuple of meshes.
    It gives access to all the vertices of all the sub-meshes as if it were a mesh itself.
    Collections can be nested to store meshes in a tree structure.

    Parameters
    ----------
    meshes: Mesh or CollectionOfMeshes

    name : str, optional
        a name for the collection
    """

    def __new__(cls, meshes, name=None):
        self = super().__new__(cls, meshes)

        for mesh in self:
            assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)

        if name is None:
            self.name = self.format_name(", ".join((mesh.name for mesh in meshes))[:-2])
        else:
            self.name = name

        return self

    def format_name(self, options_string: str) -> str:
        """Helper function to generate a name for the collection.
        Is expected to be used also in child classes."""
        if len(options_string) > NAME_MAX_LENGTH:
            options_string = options_string[:-3] + "..."
        return f"{self.__class__.__name__}({options_string})"

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def tree_view(self, **kwargs):
        body_tree_views = []
        for i, mesh in enumerate(self):
            tree_view = mesh.tree_view(**kwargs)
            if i == len(self)-1:
                prefix = ' └─'
                shift  = '   '
            else:
                prefix = ' ├─'
                shift  = ' │ '
            body_tree_views.append(prefix + tree_view.replace('\n', '\n' + shift))

        return self.name + '\n' + '\n'.join(body_tree_views)

    def copy(self, name=None):
        from copy import deepcopy
        new_mesh = deepcopy(self)
        if name is not None:
            new_mesh.name = name
        return new_mesh

    ##############
    # Properties #
    ##############

    @property
    def nb_submeshes(self):
        return len(self)

    @property
    def nb_vertices(self):
        return sum(mesh.nb_vertices for mesh in self)

    @property
    def nb_faces(self):
        return sum(mesh.nb_faces for mesh in self)

    @property
    def volume(self):
        return sum(mesh.volume for mesh in self)

    @property
    def vertices(self):
        return np.concatenate([mesh.vertices for mesh in self])

    @property
    def faces(self):
        """Return the indices of the vertices forming each of the faces. For the
        later submeshes, the indices of the vertices has to be shifted to
        correspond to their index in the concatenated array self.vertices.
        """
        nb_vertices = accumulate(chain([0], (mesh.nb_vertices for mesh in self[:-1])))
        return np.concatenate([mesh.faces + nbv for mesh, nbv in zip(self, nb_vertices)])

    @property
    def faces_normals(self):
        return np.concatenate([mesh.faces_normals for mesh in self])

    @property
    def faces_areas(self):
        return np.concatenate([mesh.faces_areas for mesh in self])

    @property
    def faces_centers(self):
        return np.concatenate([mesh.faces_centers for mesh in self])

    @property
    def faces_radiuses(self):
        return np.concatenate([mesh.faces_radiuses for mesh in self])

    def indices_of_mesh(self, mesh_index: int) -> slice:
        """Return the indices of the faces for the sub-mesh given as argument."""
        start = sum((mesh.nb_faces for mesh in self[:mesh_index]))  # Number of faces in previous meshes
        return slice(start, start + self[mesh_index].nb_faces)

    ##################
    # Transformation #
    ##################

    def merge(self, name=None) -> Mesh:
        """Merge the sub-meshes and return a full mesh.
        If the collection contains other collections, they are merged recursively.
        Optionally, a new name can be given to the resulting mesh."""
        if name is None:
            name = '+'.join(mesh.name for mesh in self)
        merged = Mesh(self.vertices, self.faces, name=name)
        merged.merge_duplicates()
        merged.heal_triangles()
        return merged

    @inplace_transformation
    def translate(self, vector):
        for mesh in self:
            mesh.translate(vector)

    @inplace_transformation
    def rotate(self, axis, angle):
        for mesh in self:
            mesh.rotate(axis, angle)

    @inplace_transformation
    def mirror(self, plane):
        for mesh in self:
            mesh.mirror(plane)

    @inplace_transformation
    def keep_immersed_part(self, **kwargs):
        for mesh in self:
            mesh.keep_immersed_part(**kwargs)
        # TODO: Prune empty meshes?

    # @inplace_transformation
    # def prune_empty_meshes(self):
    #     """Remove empty meshes from the collection."""
    #     for mesh in self:
    #         if mesh.nb_faces == 0 and mesh.nb_vertices == 0:

    def show(self):
        self.merge().show()

    def show_matplotlib(self, *args, **kwargs):
        self.merge().show_matplotlib(*args, **kwargs)
