#!/usr/bin/env python
# coding: utf-8
"""Floating bodies to be used in radiation-diffraction problems.

class FloatingBody
"""

import logging
import copy

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.geometry import Plane
from meshmagick.mesh_clipper import MeshClipper

LOG = logging.getLogger(__name__)


class FloatingBody(Mesh):
    """A floating body described as a mesh and some degrees of freedom.

    The mesh structure is inherited from meshmagick Mesh class (see
    documentation of this class for more details). The degrees of freedom
    (dofs) are stored as a dict associating a name to a 1 dimensional array of
    length equal to the number of faces in the mesh.
    """

    #######################################
    #  Initialisation and transformation  #
    #######################################

    def __init__(self, *args, **kwargs):
        Mesh.__init__(self, *args, **kwargs)
        self._compute_radiuses()
        self.nb_matrices_to_keep = 1
        self.dofs = {}
        LOG.info(f"New floating body: {self.name}.")

    def __repr__(self):
        return self.name

    @staticmethod
    def from_file(filename, file_format):
        """Create a FloatingBody from a mesh file using meshmagick."""
        from meshmagick.mmio import load_mesh
        from capytaine.symmetries import ReflectionSymmetry, xOz_Plane

        vertices, faces = load_mesh(filename, file_format)
        body = FloatingBody(vertices, faces, name=filename)

        if file_format == 'mar':
            with open(filename, 'r') as fi:
                header = fi.readline()
                _, sym = header.split()
                if int(sym) == 1:
                    body = ReflectionSymmetry(body, plane=xOz_Plane)

        return body

    def __add__(self, body_to_add):
        """Create a new CollectionOfFloatingBody from the combination of two of them."""
        from capytaine.bodies_collection import CollectionOfFloatingBodies
        return CollectionOfFloatingBodies([self, body_to_add])

    def as_FloatingBody(self, name=None):
        if name is not None:
            LOG.debug(f"Rename {self.name} as {name}.")
            self.name = name
        return self

    def copy(self, name=None):
        new_body = copy.deepcopy(self)
        if name is None:
            new_body.name = f"copy_of_{self.name}"
            LOG.debug(f"Copy {self.name}.")
        else:
            new_body.name = name
            LOG.debug(f"Copy {self.name} under the name {name}.")
        return new_body

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """Create a new FloatingBody by extracting some faces from the mesh."""
        if return_index:
            new_body, id_v = Mesh.extract_faces(self, id_faces_to_extract, return_index)
        else:
            new_body = Mesh.extract_faces(self, id_faces_to_extract, return_index)
        new_body.__class__ = FloatingBody
        new_body.nb_matrices_to_keep = self.nb_matrices_to_keep
        LOG.info(f"Extract floating body from {self.name}.")

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs[name] = dof[id_faces_to_extract]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    def get_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Remove the parts of the body above the free surface or below the sea bottom."""
        clipped_mesh = MeshClipper(self,
                                   plane=Plane(normal=(0.0, 0.0, 1.0),
                                               scalar=free_surface)).clipped_mesh

        if sea_bottom > -np.infty:
            clipped_mesh = MeshClipper(clipped_mesh,
                                       plane=Plane(normal=(0.0, 0.0, -1.0),
                                                   scalar=-sea_bottom)).clipped_mesh

        clipped_mesh.remove_unused_vertices()

        LOG.info(f"Clip floating body {self.name}.")
        return FloatingBody(clipped_mesh.vertices, clipped_mesh.faces)

    def add_translation_dof(self, direction=(1.0, 0.0, 0.0), name=None):
        if name is None:
            name = f"dof_{self.nb_dofs}_translation"
        self.dofs[name] = self.faces_normals @ direction

    def add_rotation_dof(self, axis_direction=(0.0, 0.0, 1.0), axis_point=(0.0, 0.0, 0.0), name=None):
        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        # TODO: Rewrite more efficiently and/or elegantly
        dof = np.empty((self.nb_faces, ), dtype=np.float32)
        for i, (cdg, normal) in enumerate(zip(self.faces_centers, self.faces_normals)):
            dof[i] = np.cross(axis_point - cdg, axis_direction) @ normal
        self.dofs[name] = dof

    ########################
    #  Various properties  #
    ########################

    @property
    def nb_dofs(self):
        """Number of degrees of freedom."""
        return len(self.dofs)

    @property
    def faces_radiuses(self):
        """Get the array of faces radiuses of the mesh."""
        if 'faces_radiuses' not in self.__internals__:
            self._compute_radiuses()
        return self.__internals__['faces_radiuses']

    def _compute_radiuses(self):
        """Compute the radiuses of the faces of the mesh.

        The radius is defined here as the maximal distance between the center
        of mass of a cell and one of its points."""
        from numpy.linalg import norm
        faces_radiuses = np.zeros(self.nb_faces, dtype=np.float32)
        for j in range(self.nb_faces): # TODO: optimize by array broadcasting
            faces_radiuses[j] = max(
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 0], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 1], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 2], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 3], 0:3]),
                )
        self.__internals__["faces_radiuses"] = faces_radiuses

    #######################################
    #  Computation of influence matrices  #
    #######################################

    def build_matrices(self, solver, other_body, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0, **kwargs):
        """Return the influence matrices of self on other_body."""

        LOG.debug(f"\tEvaluating matrix of {self.name} on {other_body.name} for depth={free_surface-sea_bottom:.2e} and k={wavenumber:.2e}")

        S = np.zeros((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
        V = np.zeros((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

        S0, V0 = solver._build_matrices_0(self, other_body)
        S += S0
        V += V0

        if free_surface < np.infty:

            S1, V1 = solver._build_matrices_1(self, other_body, free_surface, sea_bottom)
            S += S1
            V += V1

            S2, V2 = solver._build_matrices_2(self, other_body, free_surface, sea_bottom, wavenumber)
            S += S2
            V += V2

        return S, V
