#!/usr/bin/env python
# coding: utf-8
"""Floating bodies to be used in radiation-diffraction problems."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging
import copy
from itertools import chain, accumulate

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.mesh.meshes_collection import CollectionOfMeshes
from capytaine.tools.geometry import Abstract3DObject, inplace_or_not

from capytaine.ui.vtk.mesh_viewer import FloatingBodyViewer

LOG = logging.getLogger(__name__)

TRANSLATION_DOFS_DIRECTIONS = {"surge": (1, 0, 0), "sway": (0, 1, 0), "heave": (0, 0, 1)}
ROTATION_DOFS_AXIS = {"roll": (1, 0, 0), "pitch": (0, 1, 0), "yaw": (0, 0, 1)}


class FloatingBody(Abstract3DObject):

    def __init__(self, mesh=None, dofs=None, name=None):
        """A floating body described as a mesh and some degrees of freedom.

        The mesh structure is stored as an instance from the Mesh class of meshmagick (see
        documentation of this class for more details) or as a CollectionOfMeshes from
        capytaine.meshes_collection. The latter is a tuple of meshes or of other collections.

        The degrees of freedom (dofs) are stored as a dict associating a name to an array
        of shape (nb_faces, 3) associating a vector to each face of the body.

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
    def from_file(filename: str, file_format=None, name=None) -> 'FloatingBody':
        """Create a FloatingBody from a mesh file using meshmagick."""
        from capytaine.io.mesh_loaders import load_mesh
        if name is None:
            name = filename
        mesh = load_mesh(filename, file_format, name=f"{name}_mesh")
        return FloatingBody(mesh, name=name)

    def __lt__(self, other: 'FloatingBody') -> bool:
        """Arbitrary order. The point is to sort together the problems involving the same body."""
        return self.name < other.name

    # @property
    # def center_of_buoyancy(self):
    #     mesh = self.mesh.merge() if isinstance(self.mesh, CollectionOfMeshes) else self.mesh
    #     return Hydrostatics(mesh).buoyancy_center

    # @property
    # def displacement_volume(self):
    #     mesh = self.mesh.merge() if isinstance(self.mesh, CollectionOfMeshes) else self.mesh
    #     return Hydrostatics(mesh).displacement_volume

    @property
    def center_of_gravity(self):
        if hasattr(self, 'center'):
            return self.center

    ##########
    #  Dofs  #
    ##########

    @property
    def nb_dofs(self) -> int:
        """Number of degrees of freedom."""
        return len(self.dofs)

    def add_translation_dof(self, direction=None, name=None, amplitude=1.0):
        """Add a new translation dof (in place).
        If no direction is given, the code tries to infer it from the name.

        Parameters
        ----------
        direction : array of shape (3,), optional
            the direction of the translation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0 m/s)
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

        motion = np.empty((self.mesh.nb_faces, 3))
        motion[:, :] = direction
        self.dofs[name] = amplitude * motion

    def add_rotation_dof(self, axis=None, name=None, amplitude=1.0):
        """Add a new rotation dof (in place).
        If no axis is given, the code tries to infer it from the name.

        Parameters
        ----------
        axis: Axis, optional
            the axis of the rotation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0)
        """
        if axis is None:
            if name is not None and name.lower() in ROTATION_DOFS_AXIS:
                axis_direction = ROTATION_DOFS_AXIS[name.lower()]
                if hasattr(self, 'center_of_gravity'):
                    axis_point = self.center_of_gravity
                    LOG.info(f"The rotation dof {name} have been initialized "
                             f"around the center of gravity of {self.name}.")
                else:
                    axis_point = np.array([0, 0, 0])
                    LOG.warning(f"The rotation dof {name} have been initialized "
                                f"around the origin of the domain (0, 0, 0).")
            else:
                raise ValueError("A direction needs to be specified for the dof.")
        else:
            axis_point = axis.point
            axis_direction = axis.vector

        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        motion = np.cross(axis_point - self.mesh.faces_centers, axis_direction)
        self.dofs[name] = amplitude * motion

    def add_all_rigid_body_dofs(self):
        self.add_translation_dof(name="Surge")
        self.add_translation_dof(name="Sway")
        self.add_translation_dof(name="Heave")
        self.add_rotation_dof(name="Roll")
        self.add_rotation_dof(name="Pitch")
        self.add_rotation_dof(name="Yaw")

    ###################
    # Transformations #
    ###################

    def __add__(self, body_to_add: 'FloatingBody') -> 'FloatingBody':
        return self.join_bodies(body_to_add)

    def join_bodies(*bodies, name=None) -> 'FloatingBody':
        if name is None:
            name = "+".join(body.name for body in bodies)
        meshes = CollectionOfMeshes([body.mesh for body in bodies], name=f"{name}_mesh")
        dofs = FloatingBody.combine_dofs(bodies)
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
                new_dof = np.zeros((total_nb_faces, 3))
                new_dof[nbf:nbf+len(dof), :] = dof
                dofs['_'.join([body.name, name])] = new_dof
        return dofs

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
            new_body.dofs[name] = dof[id_faces_to_extract, :]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    def get_immersed_part(self, name=None, **kwargs):
        """Return a body for which the parts of the mesh above the free surface or below the sea
        bottom have been removed.
        Dofs are lost in the process.
        TODO: Also clip dofs.
        """
        new_body_mesh = self.mesh.get_immersed_part(**kwargs)
        if new_body_mesh is None:
            raise Exception(f"Trying to clip the mesh of {self.name}, but it does not have any wet faces")

        if name is None:
            name = f"{self.name}_clipped"
        LOG.info(f"Clip floating body {self.name} to create {name}.")

        return FloatingBody(new_body_mesh, name=name)

    @inplace_or_not
    def mirror(self, *args):
        self.mesh.mirror(*args)
        # TODO: Also mirror dofs
        LOG.warning(f"The dofs of {self.name} have not been mirrored with its mesh.")
        return

    @inplace_or_not
    def translate(self, *args):
        if hasattr(self, 'center'):
            self.center += args[0]
        self.mesh.translate(*args)
        return

    @inplace_or_not
    def rotate(self, *args):
        # TODO: Also rotate dofs
        LOG.warning(f"The dofs of {self.name} have not been rotated with its mesh.")
        return self.mesh.rotate(*args)

    #############
    #  Display  #
    #############

    def __str__(self):
        return self.name

    def __repr__(self):
        return (f"{self.__class__.__name__}(mesh={self.mesh.name}, "
                f"dofs={{{', '.join(self.dofs.keys())}}}, name={self.name})")

    def show(self, dof=None):
        # TODO: Broken. Rewrite with display of dofs.
        import vtk
        vtk_polydata = self.mesh._vtk_polydata()

        if dof is not None:
            vtk_data_array = vtk.vtkFloatArray()
            vtk_data_array.SetNumberOfComponents(3)
            vtk_data_array.SetNumberOfTuples(self.mesh.nb_faces)
            for i, vector in enumerate(self.dofs[dof]):
                vtk_data_array.SetTuple3(i, *vector)
            vtk_polydata.GetCellData().SetVectors(vtk_data_array)

        viewer = FloatingBodyViewer()
        viewer.add_polydata(vtk_polydata)
        viewer.show()
        viewer.finalize()

    def show_matplotlib(self, *args, **kwargs):
        return self.mesh.show_matplotlib(*args, **kwargs)

