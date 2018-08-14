#!/usr/bin/env python
# coding: utf-8
"""Storing a set of meshes."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import copy
from itertools import chain, accumulate

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.tools.geometry import Plane

NAME_MAX_LENGTH = 180


class CollectionOfMeshes:
    name: str
    # submeshes: Iterable[Union[Mesh, CollectionOfMeshes]]

    def __init__(self, meshes, name=None):
        """A list of meshes.
        It gives access to all the vertices of all the sub-meshes as if it were a mesh itself.
        Collections can be nested to store meshes in a tree structure.

        Parameters
        ----------
        meshes : iterable of Mesh or of CollectionOfMeshes
            the meshes contained in the collection
        name : str
            a name for the collection (optional)
        """
        self.submeshes = tuple(meshes)
        for mesh in self.submeshes:
            assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)

        if name is None:
            self.name = self.format_name(", ".join((mesh.name for mesh in meshes))[:-2])
        else:
            self.name = name

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

    def copy(self, name=None):
        new_mesh = copy.deepcopy(self)
        if name is not None:
            new_mesh.name = name
        return new_mesh

    ##############
    # Properties #
    ##############

    @property
    def nb_submeshes(self):
        return len(self.submeshes)

    @property
    def nb_vertices(self):
        return sum(mesh.nb_vertices for mesh in self.submeshes)

    @property
    def nb_faces(self):
        return sum(mesh.nb_faces for mesh in self.submeshes)

    @property
    def volume(self):
        return sum(mesh.volume for mesh in self.submeshes)

    @property
    def vertices(self):
        return np.concatenate([mesh.vertices for mesh in self.submeshes])

    @property
    def faces(self):
        """Return the indices of the vertices forming each of the faces. For the
        later subbodies, the indices of the vertices has to be shifted to
        correspond to their index in the concatenated array self.vertices.
        """
        nb_vertices = accumulate(chain([0], (mesh.nb_vertices for mesh in self.submeshes[:-1])))
        return np.concatenate([mesh.faces + nbv for mesh, nbv in zip(self.submeshes, nb_vertices)])

    @property
    def faces_normals(self):
        return np.concatenate([mesh.faces_normals for mesh in self.submeshes])

    @property
    def faces_areas(self):
        return np.concatenate([mesh.faces_areas for mesh in self.submeshes])

    @property
    def faces_centers(self):
        return np.concatenate([mesh.faces_centers for mesh in self.submeshes])

    @property
    def faces_radiuses(self):
        return np.concatenate([mesh.faces_radiuses for mesh in self.submeshes])

    def indices_of_mesh(self, mesh_index: int) -> slice:
        """Return the indices of the faces for the sub-mesh given as argument."""
        start = sum((mesh.nb_faces for mesh in self.submeshes[:mesh_index]))  # Number of faces in previous meshes
        return slice(start, start + self.submeshes[mesh_index].nb_faces)

    ##################
    # Transformation #
    ##################

    def merge(self, name: str=None) -> Mesh:
        """Merge the sub-meshes and return a full mesh.
        If the collection contains other collections, they are merged recursively.
        Optionally, a new name can be given to the resulting mesh."""
        components = (mesh if isinstance(mesh, Mesh) else mesh.merge() for mesh in self.submeshes)
        init = next(components)
        merged = sum(components, init)
        merged.merge_duplicates()
        merged.heal_triangles()
        if name is not None:
            merged.name = name
        return merged

    def get_immersed_part(self, **kwargs):
        clipped_meshes = []
        for mesh in self.submeshes:
            m = mesh.get_immersed_part(**kwargs)
            if m is not None:
                clipped_meshes.append(m)

        if len(clipped_meshes) > 0:
            return CollectionOfMeshes(clipped_meshes)
        else:
            return None

    def mirror(self, plane):
        for mesh in self.submeshes:
            mesh.mirror(plane)

    def translate_x(self, value):
        for mesh in self.submeshes:
            mesh.translate_x(value)

    def translate_y(self, value):
        for mesh in self.submeshes:
            mesh.translate_y(value)

    def translate_z(self, value):
        for mesh in self.submeshes:
            mesh.translate_z(value)

    def translate(self, vector):
        for mesh in self.submeshes:
            mesh.translate(vector)

    def rotate_x(self, value):
        for mesh in self.submeshes:
            mesh.rotate_x(value)

    def rotate_y(self, value):
        for mesh in self.submeshes:
            mesh.rotate_y(value)

    def rotate_z(self, value):
        for mesh in self.submeshes:
            mesh.rotate_z(value)

    def rotate(self, vector):
        for mesh in self.submeshes:
            mesh.rotate(vector)

    def _vtk_polydata(self):
        return self.merge()._vtk_polydata()

    def show(self):
        self.merge().show()

    def show_matplotlib(self, *args, **kwargs):
        self.merge().show_matplotlib(*args, **kwargs)
