#!/usr/bin/env python
# coding: utf-8

import logging
from itertools import chain, accumulate

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import FloatingBody


LOG = logging.getLogger(__name__)


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

        # Name of the body collection.
        self.name = "union_of_{name_list}_and_{last_body_name}".format(
                name_list='_'.join((body.name for body in bodies[:-1])),
                last_body_name=bodies[-1].name
                )
        LOG.info(f"New body: {self.name}.")

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

    def as_FloatingBody(self):
        """Merge the mesh of the bodies of the collection into one mesh."""
        new_body = self.subbodies[0].as_FloatingBody().copy()
        for body in self.subbodies[1:]:
            new_body = Mesh.__add__(new_body, body.as_FloatingBody())
        new_body.merge_duplicates()
        new_body.heal_triangles()
        new_body.__class__ = FloatingBody
        new_body.name = self.name
        new_body.dofs = self.dofs
        new_body.nb_matrices_to_keep = 1
        return new_body

    def __add__(self, body_to_add):
        if isinstance(body_to_add, CollectionOfFloatingBodies):
            return CollectionOfFloatingBodies(self.subbodies + body_to_add.subbodies)
        else:
            return CollectionOfFloatingBodies(self.subbodies + [body_to_add])

    ########################
    #  Various properties  #
    ########################

    def __str__(self):
        return self.name

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
        start = sum((body.nb_faces for body in self.subbodies[:body_index]))
        return slice(start, start + self.subbodies[body_index].nb_faces)

    #######################################
    #  Computation of influence matrices  #
    #######################################

    def build_matrices(self, other_body, **kwargs):
        """Return the influence matrices of self on other body."""
        LOG.debug(f"Evaluating matrix of {self.name} on {other_body.name}.")

        S = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
        V = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

        nb_faces = list(accumulate(chain([0], (body.nb_faces for body in self.subbodies))))
        for (i, j), body in zip(zip(nb_faces, nb_faces[1:]), self.subbodies):
            matrix_slice = (slice(i, j), slice(None, None))
            S[matrix_slice], V[matrix_slice] = body.build_matrices(other_body, **kwargs)

        return S, V

