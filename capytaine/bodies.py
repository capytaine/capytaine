#!/usr/bin/env python
# coding: utf-8
"""Floating bodies to be used in radiation-diffraction problems.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
import copy

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.geometry import Plane
from meshmagick.mesh_clipper import MeshClipper

from capytaine.meshes_collection import CollectionOfMeshes

LOG = logging.getLogger(__name__)


class FloatingBody:
    """A floating body described as a mesh and some degrees of freedom.

    The mesh structure is stored in a meshmagick Mesh (see documentation of
    this class for more details). The degrees of freedom (dofs) are stored as a
    dict associating a name to a 1 dimensional array of length equal to the
    number of faces in the mesh.
    """

    def __init__(self, mesh=None, name=None):
        if mesh is None:
            self.mesh = Mesh(np.zeros((0, 3)), np.zeros((0, 4)))  # Dummy mesh
        else:
            assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)
            self.mesh = mesh

        if name is None:
            self.name = self.mesh.name
        else:
            self.name = name

        self.dofs = {}

        LOG.info(f"New floating body: {self.name}.")

    @staticmethod
    def from_file(filename, file_format="mar"):
        """Create a FloatingBody from a mesh file using meshmagick."""
        from meshmagick.mmio import load_mesh
        from meshmagick.geometry import xOz_Plane
        from capytaine.symmetries import ReflectionSymmetry

        vertices, faces = load_mesh(filename, file_format)
        body = FloatingBody(Mesh(vertices, faces, name=f"{filename}_mesh"), name=filename)

        if file_format == 'mar':
            with open(filename, 'r') as fi:
                header = fi.readline()
                _, sym = header.split()
                if int(sym) == 1:
                    body = ReflectionSymmetry(body, plane=xOz_Plane)

        return body

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __lt__(self, other):
        """Arbitrary order. The point is to sort together the problems involving the same body."""
        return self.name < other.name

    # def __add__(self, body_to_add):
    #     """Create a new CollectionOfFloatingBody from the combination of two FloatingBodies."""
    #     return FloatingBody(CollectionOfMeshes([self.mesh, self.CollectionOfFloatingBodies([self, body_to_add])

    ##########
    #  Dofs  #
    ##########

    @property
    def nb_dofs(self):
        """Number of degrees of freedom."""
        return len(self.dofs)

    def add_translation_dof(self, direction=(1.0, 0.0, 0.0), name=None):
        """Helper to define a new translation dof."""
        if name is None:
            name = f"dof_{self.nb_dofs}_translation"
        self.dofs[name] = self.mesh.faces_normals @ direction

    def add_rotation_dof(self, axis_direction=(0.0, 0.0, 1.0), axis_point=(0.0, 0.0, 0.0), name=None):
        """Helper to define a new rotation dof."""
        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        axis_direction = np.asarray(axis_direction)
        axis_point = np.asarray(axis_point)

        motion = np.cross(axis_point - self.mesh.faces_centers, axis_direction)
        motion[motion != 0] /= np.linalg.norm(motion[motion != 0])
        dof = np.sum(motion * self.mesh.faces_normals, axis=1)

        self.dofs[name] = dof

    ###################
    # Transformations #
    ###################

    def copy(self, name=None):
        new_body = copy.deepcopy(self)
        if name is None:
            new_body.name = f"copy_of_{self.name}"
            LOG.debug(f"Copy {self.name}.")
        else:
            new_body.name = name
            LOG.debug(f"Copy {self.name} under the name {name}.")
        return new_body

    def get_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Remove the parts of the body above the free surface or below the sea bottom.
        Dofs are lost in the process."""
        # TODO: Also clip dofs
        clipped_mesh = MeshClipper(self.mesh,
                                   plane=Plane(normal=(0.0, 0.0, 1.0),
                                               scalar=free_surface)).clipped_mesh

        if sea_bottom > -np.infty:
            clipped_mesh = MeshClipper(clipped_mesh,
                                       plane=Plane(normal=(0.0, 0.0, -1.0),
                                                   scalar=-sea_bottom)).clipped_mesh

        clipped_mesh.remove_unused_vertices()
        LOG.info(f"Clip floating body {self.name}.")
        return FloatingBody(clipped_mesh, name=f"{self.name}_clipped")

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """Create a new FloatingBody by extracting some faces from the mesh."""
        # TODO: CollectionOfMeshes
        if return_index:
            new_mesh, id_v = Mesh.extract_faces(self.mesh, id_faces_to_extract, return_index)
        else:
            new_mesh = Mesh.extract_faces(self.mesh, id_faces_to_extract, return_index)
        new_body = FloatingBody(new_mesh)
        LOG.info(f"Extract floating body from {self.name}.")

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs[name] = dof[id_faces_to_extract]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    def mirror(self, *args):
        return self.mesh.mirror(*args)

    def translate_x(self, *args):
        return self.mesh.translate_x(*args)

    def translate_y(self, *args):
        return self.mesh.translate_y(*args)

    def translate_z(self, *args):
        return self.mesh.translate_z(*args)

    def translate(self, *args):
        return self.mesh.translate(*args)

    def rotate_x(self, *args):
        # TODO: Also rotate dofs
        return self.mesh.rotate_x(*args)

    def rotate_y(self, *args):
        # TODO: Also rotate dofs
        return self.mesh.rotate_y(*args)

    def rotate_z(self, *args):
        # TODO: Also rotate dofs
        return self.mesh.rotate_z(*args)

    def rotate(self, *args):
        # TODO: Also rotate dofs
        return self.mesh.rotate(*args)

    def show(self):
        # TODO: Also show dofs
        return self.mesh.show()

    def show_matplotlib(self, *args, **kwargs):
        return self.mesh.show_matplotlib(*args, **kwargs)

#     @staticmethod
#     def repeat_dof(bodies):
#         """Combine the degrees of freedom of the subbodies."""
#         dofs = {}
#         cum_nb_faces = accumulate(chain([0], (body.mesh.nb_faces for body in bodies)))
#         total_nb_faces = sum(body.mesh.nb_faces for body in bodies)
#         for body, nbf in zip(bodies, cum_nb_faces):
#             # nbf is the cumulative number of faces of the previous subbodies,
#             # that is the offset of the indices of the faces of the current body.
#             for name, dof in body.dofs.items():
#                 dofs['_'.join([body.name, name])] = np.r_[np.zeros(nbf),
#                                                           dof,
#                                                           np.zeros(total_nb_faces - len(dof) - nbf)]
#         return dofs
