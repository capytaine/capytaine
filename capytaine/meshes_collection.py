#!/usr/bin/env python
# coding: utf-8

import copy
from typing import TypeVar

from itertools import chain, accumulate

import numpy as np

from meshmagick.mesh import Mesh

NAME_MAX_LENGTH = 180


class CollectionOfMeshes:
    """A list of meshes"""

    def __init__(self, meshes, name=None):
        for mesh in meshes:
            assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)
        self.submeshes = meshes

        if name is None:
            self.name = self.format_name(", ".join((mesh.name for mesh in meshes))[:-2])
        else:
            self.name = name

    def format_name(self, options_string):
        if len(options_string) > NAME_MAX_LENGTH:
            options_string = options_string[:-3] + "..."
        return f"{self.__class__.__name__}({options_string})"

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def copy(self):
        return copy.deepcopy(self)

    # Properties

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

    def indices_of_mesh(self, mesh_index):
        start = sum((mesh.nb_faces for mesh in self.submeshes[:mesh_index]))
        return slice(start, start + self.submeshes[mesh_index].nb_faces)

    # Transformation

    def merge(self) -> Mesh:
        components = (mesh.merge() if isinstance(mesh, CollectionOfMeshes) else mesh for mesh in self.submeshes)  # Ensure components have been merged
        init = next(components)
        merged = sum(components, init)
        merged.merge_duplicates()
        merged.heal_triangles()
        return merged

    def mirror(self, plane):
        for mesh in self.submeshes:
            mesh.mirror(plane)
        return

    def translate_x(self, value):
        for mesh in self.submeshes:
            mesh.translate_x(value)
        return

    def translate_y(self, value):
        for mesh in self.submeshes:
            mesh.translate_y(value)
        return

    def translate_z(self, value):
        for mesh in self.submeshes:
            mesh.translate_z(value)
        return

    def translate(self, vector):
        for mesh in self.submeshes:
            mesh.translate(vector)
        return

    def rotate_x(self, value):
        for mesh in self.submeshes:
            mesh.rotate_x(value)
        return

    def rotate_y(self, value):
        for mesh in self.submeshes:
            mesh.rotate_y(value)
        return

    def rotate_z(self, value):
        for mesh in self.submeshes:
            mesh.rotate_z(value)
        return

    def rotate(self, vector):
        for mesh in self.submeshes:
            mesh.rotate(vector)
        return

    def show(self):
        self.merge().show()


MeshType = TypeVar('MeshType', Mesh, CollectionOfMeshes)
