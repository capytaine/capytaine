#!/usr/bin/env python
# coding: utf-8

import logging
from itertools import chain, accumulate

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import CMesh, FloatingBody


LOG = logging.getLogger(__name__)

NAME_MAX_LENGTH = 180


class CollectionOfFloatingBodies(FloatingBody):
    """A body composed of several floating bodies."""

    #######################################
    #  Initialisation and transformation  #
    #######################################

    def __init__(self, bodies):
        """Initialize the body."""

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

        new_body = self.subbodies[0].as_FloatingBody().copy(name=name)
        for body in self.subbodies[1:]:
            new_body.mesh = Mesh.__add__(new_body.mesh, body.as_FloatingBody().mesh)
            LOG.debug(f"Add mesh of {body.name} to {name}.")
        new_body.mesh.__class__ = CMesh
        new_body.mesh.merge_duplicates()
        new_body.mesh.heal_triangles()
        new_body.name = name
        new_body.dofs = self.dofs  # TODO: is broken if the subbodies have faces in common.
        new_body.nb_matrices_to_keep = 1
        new_body.__internals__ = {}
        LOG.info(f"Merged collection of bodies {self.name} into floating body {new_body.name}.")
        return new_body

    def __add__(self, body_to_add):
        if isinstance(body_to_add, CollectionOfFloatingBodies):
            return CollectionOfFloatingBodies(self.subbodies + body_to_add.subbodies)
        else:
            return CollectionOfFloatingBodies(self.subbodies + [body_to_add])

    ########################
    #  Various properties  #
    ########################

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
        """Return the indices of the verices forming each of the faces. For the
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

    def indices_of_body(self, body_index):
        start = sum((body.mesh.nb_faces for body in self.subbodies[:body_index]))
        return slice(start, start + self.subbodies[body_index].mesh.nb_faces)

    #######################################
    #  Computation of influence matrices  #
    #######################################

    def build_matrices(self, solver, other_body, **kwargs):
        """Return the influence matrices of self on other body."""
        LOG.debug(f"Evaluating matrix of {self.name} on {other_body.name}.")

        S = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
        V = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

        nb_faces = list(accumulate(chain([0], (body.nb_faces for body in self.subbodies))))
        for (i, j), body in zip(zip(nb_faces, nb_faces[1:]), self.subbodies):
            matrix_slice = (slice(i, j), slice(None, None))
            S[matrix_slice], V[matrix_slice] = body.build_matrices(solver, other_body, **kwargs)

        return S, V

