#!/usr/bin/env python
# coding: utf-8

import logging
from itertools import chain, accumulate

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import FloatingBody


LOG = logging.getLogger(__name__)

NAME_MAX_LENGTH = 180


class CollectionOfFloatingBodies(FloatingBody):
    """A body composed of several floating bodies."""

    #######################################
    #  Initialisation and transformation  #
    #######################################

    def __init__(self, bodies):
        """Initialize the body.

        Parameters
        ----------
        bodies: list (or any iterable) of FloatingBody
            the bodies composing the collection
        """

        for body in bodies:
            assert isinstance(body, FloatingBody)

        self.subbodies = bodies

        names_of_subbodies = ', '.join(body.name for body in self.subbodies)
        if len(names_of_subbodies) > NAME_MAX_LENGTH:
            names_of_subbodies = names_of_subbodies[:NAME_MAX_LENGTH-3] + "..."
        self.name = f"CollectionOfFloatingBodies([{names_of_subbodies}])"

        LOG.debug(f"New collection of bodies: {self.name}.")

        # Combine the degrees of freedom of the subbodies.
        self.dofs = {}
        cum_nb_faces = accumulate(chain([0], (body.nb_faces for body in self.subbodies)))
        total_nb_faces = sum(body.nb_faces for body in self.subbodies)
        for body, nbf in zip(bodies, cum_nb_faces):
            # nbf is the cumulative number of faces of the previous subbodies,
            # that is the offset of the indices of the faces of the current body.
            for name, dof in body.dofs.items():
                self.dofs['_'.join([body.name, name])] = np.r_[
                        np.zeros(nbf),
                        dof,
                        np.zeros(total_nb_faces - len(dof) - nbf),
                        ]

    def as_FloatingBody(self, name=None):
        """Merge the mesh of the bodies of the collection into one mesh."""

        if name is None:
            name = f"union_of_{'_and_'.join((body.name for body in self.subbodies))}"
            if len(name) > NAME_MAX_LENGTH:
                name = name[:NAME_MAX_LENGTH-3] + "..."

        new_mesh = self.subbodies[0].as_FloatingBody().mesh
        for body in self.subbodies[1:]:
            new_mesh = Mesh.__add__(new_mesh, body.as_FloatingBody().mesh)
            LOG.debug(f"Add mesh of {body.name} to {name}.")
        new_mesh.merge_duplicates()
        new_mesh.heal_triangles()
        new_body = FloatingBody(new_mesh, name=name)
        new_body.dofs = self.dofs  # TODO: is broken if the subbodies have faces in common.
        LOG.info(f"Merged collection of bodies {self.name} into floating body {new_body.name}.")
        return new_body

    def __add__(self, body_to_add):
        if isinstance(body_to_add, CollectionOfFloatingBodies):
            return CollectionOfFloatingBodies(self.subbodies + body_to_add.subbodies)
        else:
            return CollectionOfFloatingBodies(self.subbodies + [body_to_add])

    @property
    def nb_matrices_to_keep(self):
        return max([body.nb_matrices_to_keep for body in self.subbodies])

    @nb_matrices_to_keep.setter
    def nb_matrices_to_keep(self, value):
        for body in self.subbodies:
            body.nb_matrices_to_keep = value

    @property
    def nb_subbodies(self):
        return len(self.subbodies)

    #######################
    #  Interface to mesh  #
    #######################

    @property
    def mesh(self):
        meshes = [body.mesh for body in self.subbodies]
        return sum(meshes[1:], meshes[0])

    @property
    def nb_vertices(self):
        return sum(body.nb_vertices for body in self.subbodies)

    @property
    def nb_faces(self):
        return sum(body.nb_faces for body in self.subbodies)

    @property
    def volume(self):
        return sum(body.volume for body in self.subbodies)

    @property
    def vertices(self):
        return np.concatenate([body.vertices for body in self.subbodies])

    @property
    def faces(self):
        """Return the indices of the vertices forming each of the faces. For the
        later subbodies, the indices of the vertices has to be shifted to
        correspond to their index in the concatenated array self.vertices.
        """
        nb_vertices = accumulate(chain([0], (body.nb_vertices for body in self.subbodies[:-1])))
        return np.concatenate([body.faces + nbv for body, nbv in zip(self.subbodies, nb_vertices)])

    @property
    def faces_normals(self):
        return np.concatenate([body.faces_normals for body in self.subbodies])

    @property
    def faces_areas(self):
        return np.concatenate([body.faces_areas for body in self.subbodies])

    @property
    def faces_centers(self):
        return np.concatenate([body.faces_centers for body in self.subbodies])

    @property
    def faces_radiuses(self):
        return np.concatenate([body.faces_radiuses for body in self.subbodies])

    def mirror(self, plane):
        for body in self.subbodies:
            body.mirror(plane)
        return

    def translate_x(self, value):
        for body in self.subbodies:
            body.translate_x(value)
        return

    def translate_y(self, value):
        for body in self.subbodies:
            body.translate_y(value)
        return

    def translate_z(self, value):
        for body in self.subbodies:
            body.translate_z(value)
        return

    def translate(self, vector):
        for body in self.subbodies:
            body.translate(vector)
        return

    def rotate_x(self, value):
        for body in self.subbodies:
            body.rotate_x(value)
        return

    def rotate_y(self, value):
        for body in self.subbodies:
            body.rotate_y(value)
        return

    def rotate_z(self, value):
        for body in self.subbodies:
            body.rotate_z(value)
        return

    def rotate(self, vector):
        for body in self.subbodies:
            body.rotate(vector)
        return

    def get_immersed_part(self, **kwargs):
        return CollectionOfFloatingBodies([body.get_immersed_part(**kwargs) for body in self.subbodies])

    def indices_of_body(self, body_index):
        start = sum((body.mesh.nb_faces for body in self.subbodies[:body_index]))
        return slice(start, start + self.subbodies[body_index].mesh.nb_faces)

