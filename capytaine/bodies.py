#!/usr/bin/env python
# coding: utf-8
"""Floating bodies to be used in radiation-diffraction problems.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
import copy
from itertools import chain, accumulate

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.mmio import load_mesh
from meshmagick.mesh_clipper import MeshClipper
from meshmagick.geometry import xOz_Plane, Plane

from capytaine.meshes_collection import CollectionOfMeshes
from capytaine.symmetries import ReflectionSymmetry

LOG = logging.getLogger(__name__)

TRANSLATION_DOFS_DIRECTIONS = {"surge": (1, 0, 0), "sway": (0, 1, 0), "heave": (0, 0, 1)}
ROTATION_DOFS_AXIS = {"roll": (1, 0, 0), "pitch": (0, 1, 0), "yaw": (0, 0, 1)}


class FloatingBody:

    def __init__(self, mesh=None, dofs=None, name=None):
        """A floating body described as a mesh and some degrees of freedom.

        The mesh structure is stored as an instance from the Mesh class of meshmagick (see
        documentation of this class for more details) or as a CollectionOfMeshes from
        capytaine.meshes_collection. The latter is a tuple of meshes or of other collections.

        The degrees of freedom (dofs) are stored as a dict associating a name to a 1 dimensional array
        of length equal to the number of faces in the mesh.

        Parameters
        ----------
        mesh : Mesh or CollectionOfMeshes, optional
            the mesh describing the geometry of the floating body.
            If none is given, a empty one is created.
        dofs : dict, optional
            the degrees of freedom of the body.
            If none is given, a empty dictionary is initialized.
        name : str, optional
            a name for the body.
            If none is given, the one of the mesh is used.
        """
        if mesh is None:
            mesh = Mesh(np.zeros((0, 3)), np.zeros((0, 4)), name="dummy")

        if dofs is None:
            dofs = {}

        if name is None:
            name = mesh.name

        assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)
        self.mesh = mesh
        self.dofs = dofs
        self.name = name

        LOG.info(f"New floating body: {self.name}.")

    @staticmethod
    def from_file(filename: str, file_format: str ="mar") -> 'FloatingBody':
        """Create a FloatingBody from a mesh file using meshmagick."""

        vertices, faces = load_mesh(filename, file_format)
        mesh = Mesh(vertices, faces, name=f"{filename}_mesh")

        if file_format == "mar":
            with open(filename, 'r') as fi:
                header = fi.readline()
                _, sym = header.split()
                if int(sym) == 1:
                    mesh = ReflectionSymmetry(mesh, plane=xOz_Plane)

        return FloatingBody(mesh, name=filename)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __lt__(self, other: 'FloatingBody') -> bool:
        """Arbitrary order. The point is to sort together the problems involving the same body."""
        return self.name < other.name

    def __add__(self, body_to_add: 'FloatingBody') -> 'FloatingBody':
        """Create a new CollectionOfFloatingBody from the combination of two FloatingBodies."""
        return FloatingBody.join_bodies([self, body_to_add])

    @staticmethod
    def join_bodies(bodies) -> 'FloatingBody':
        meshes = CollectionOfMeshes([body.mesh for body in bodies])
        dofs = FloatingBody.combine_dofs(bodies)
        name = name="+".join(body.name for body in bodies)
        return FloatingBody(mesh=meshes, dofs=dofs, name=name)

    @staticmethod
    def combine_dofs(bodies) -> dict:
        """Combine the degrees of freedom of several bodies."""
        dofs = {}
        cum_nb_faces = accumulate(chain([0], (body.mesh.nb_faces for body in bodies)))
        total_nb_faces = sum(body.mesh.nb_faces for body in bodies)
        for body, nbf in zip(bodies, cum_nb_faces):
            # nbf is the cumulative number of faces of the previous subbodies,
            # that is the offset of the indices of the faces of the current body.
            for name, dof in body.dofs.items():
                dofs['_'.join([body.name, name])] = np.r_[np.zeros(nbf), dof,
                                                          np.zeros(total_nb_faces - len(dof) - nbf)]
        return dofs

    ##########
    #  Dofs  #
    ##########

    @property
    def nb_dofs(self) -> int:
        """Number of degrees of freedom."""
        return len(self.dofs)

    def add_translation_dof(self, direction=None, name=None):
        """Add a new translation dof (in place).
        If no direction is given, the code tries to infer it from the name.

        Parameters
        ----------
        direction : array of shape (3,), optional
            the direction of the translation
        name : str, optional
            a name for the degree of freedom
        """
        if direction is None:
            if name is not None and name.lower() in TRANSLATION_DOFS_DIRECTIONS:
                direction = TRANSLATION_DOFS_DIRECTIONS[name.lower()]
            else:
                raise ValueError("A direction needs to be specified for the dof.")

        if name is None:
            name = f"dof_{self.nb_dofs}_translation"

        direction = np.asarray(direction)
        assert direction.shape == (3,)

        self.dofs[name] = self.mesh.faces_normals @ direction

    def add_rotation_dof(self, axis_point=np.array((0.0, 0.0, 0.0)),
                         axis_direction=None, name=None):
        """Add a new rotation dof (in place).
        If no axis direction is given, the code tries to infer it from the name.

        Parameters
        ----------
        axis_point : array of shape (3,)
            a point on the rotation axis
        axis_direction : array of shape (3,), optional
            vector directing the rotation axis
            TODO: Use Axis class?
        name : str, optional
            a name for the degree of freedom
        """
        if axis_direction is None:
            if name is not None and name.lower() in ROTATION_DOFS_AXIS:
                axis_direction = ROTATION_DOFS_AXIS[name.lower()]
            else:
                raise ValueError("A direction needs to be specified for the dof.")
        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        axis_direction = np.asarray(axis_direction)
        axis_point = np.asarray(axis_point)

        assert axis_direction.shape == (3,)
        assert axis_point.shape == (3,)

        motion = np.cross(axis_point - self.mesh.faces_centers, axis_direction)
        motion[motion != 0] /= np.linalg.norm(motion[motion != 0])  # Normalize
        dof = np.sum(motion * self.mesh.faces_normals, axis=1)

        self.dofs[name] = dof

    ###################
    # Transformations #
    ###################

    def copy(self, name=None) -> 'FloatingBody':
        """Return a deep copy of the body.

        Parameters
        ----------
        name : str, optional
            a name for the new copy
        """
        new_body = copy.deepcopy(self)
        if name is None:
            new_body.name = f"copy_of_{self.name}"
            LOG.debug(f"Copy {self.name}.")
        else:
            new_body.name = name
            LOG.debug(f"Copy {self.name} under the name {name}.")
        return new_body

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """Create a new FloatingBody by extracting some faces from the mesh.
        The dofs evolve accordingly.
        """
        if isinstance(self.mesh, CollectionOfMeshes):
            raise NotImplemented()  # TODO

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

    def get_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Return a body for which the parts of the mesh above the free surface or below the sea
        bottom have been removed.
        Dofs are lost in the process.
        TODO: Also clip dofs.
        """
        if isinstance(self.mesh, CollectionOfMeshes):
            raise NotImplemented()  # TODO

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

