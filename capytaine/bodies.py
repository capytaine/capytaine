#!/usr/bin/env python
# coding: utf-8
"""Floatting bodies to be used in radiation-diffraction problems.

class FloattingBody
"""

import numpy as np

from meshmagick.mesh import Mesh

import capytaine._Green as _Green


class FloattingBody(Mesh):
    """A floatting body composed of a mesh (inherited from Meshmagick) and
    several degrees of freedom (dof)."""

    def __init__(self, *args, **kwargs):
        Mesh.__init__(self, *args, **kwargs)
        self.dofs = {}

    @staticmethod
    def from_file(filename, file_format):
        from meshmagick.mmio import load_mesh

        vertices, faces = load_mesh(filename, file_format)

        return FloattingBody(vertices, faces, name="mesh_from_"+filename)

    def __add__(self, body_to_add):
        new_body = Mesh.__add__(self, body_to_add)
        new_body.__class__ = FloattingBody

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs['_'.join([name, self.name])] = np.r_[dof, np.zeros(body_to_add.nb_faces)]
        for name, dof in body_to_add.dofs.items():
            new_body.dofs['_'.join([name, body_to_add.name])] = np.r_[np.zeros(self.nb_faces), dof]
        return new_body

    def extract_faces(self, id_faces_to_extract, return_index=False):
        if return_index:
            new_body, id_v = Mesh.extract_faces(self, id_faces_to_extract, return_index)
        else:
            new_body = Mesh.extract_faces(self, id_faces_to_extract, return_index)
        new_body.__class__ = FloattingBody

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs[name] = dof[id_faces_to_extract]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    @property
    def nb_dofs(self):
        return len(self.dofs)

    @property
    def faces_radiuses(self):
        """Get the array of faces radiuses of the mesh

        Returns
        -------
        ndarray
        """
        if 'faces_radiuses' not in self.__internals__:
            self._compute_radiuses()
        return self.__internals__['faces_radiuses']

    def _compute_radiuses(self):
        """Update face radiuses"""
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

    def _build_matrices_0(self, body):
        """Compute the influence of self on body"""

        if 'Green0' not in self.__internals__:
            self.__internals__['Green0'] = {}

        if body not in self.__internals__['Green0']:
            S0, V0 = _Green.green_1.build_matrix_0(
                self.faces_centers, self.faces_normals,
                body.vertices,      body.faces + 1,
                body.faces_centers, body.faces_normals,
                body.faces_areas,   body.faces_radiuses,
                )
            self.__internals__['Green0'][body] = (S0, V0)

        return self.__internals__['Green0'][body]

    def _build_matrices_1(self, body, free_surface, sea_bottom):
        if 'Green1' not in self.__internals__:
            self.__internals__['Green1'] = {}
        if body not in self.__internals__['Green1']:
            self.__internals__['Green1'][body] = {}

        depth = free_surface-sea_bottom
        if depth not in self.__internals__['Green1'][body]:
            def reflect_vector(x):
                y = x.copy()
                y[:, 2] = -x[:, 2]
                return y

            if depth == np.infty:
                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*free_surface - x[:, 2] 
                    return y
            else:
                def reflect_point(x):
                    y = x.copy()
                    y[:, 2] = 2*sea_bottom - x[:, 2]
                    return y

            S1, V1 = _Green.green_1.build_matrix_0(
                reflect_point(self.faces_centers), reflect_vector(self.faces_normals),
                body.vertices,      body.faces + 1,
                body.faces_centers, body.faces_normals,
                body.faces_areas,   body.faces_radiuses,
                )

            if depth == np.infty:
                self.__internals__['Green1'][body][np.infty] = (-S1, -V1)
            else:
                self.__internals__['Green1'][body][depth] = (S1, V1)

        return self.__internals__['Green1'][body][depth]

    def _build_matrices_2(self, body, free_surface, sea_bottom, wavenumber):
        if 'Green2' not in self.__internals__:
            self.__internals__['Green2'] = {}
        if body not in self.__internals__['Green2']:
            self.__internals__['Green2'][body] = {}

        depth = free_surface - sea_bottom
        if (depth, wavenumber) not in self.__internals__['Green2'][body]:
            if depth == np.infty:
                S2, V2 = _Green.green_2.build_matrix_2(
                    self.faces_centers, self.faces_normals,
                    body.faces_centers, body.faces_areas,  
                    wavenumber,       0.0
                    )
            else:
                S2, V2 = _Green.green_2.build_matrix_2(
                    self.faces_centers, self.faces_normals,
                    body.faces_centers, body.faces_areas,  
                    wavenumber,         depth
                    )
            self.__internals__['Green2'][body][(depth, wavenumber)] = (S2, V2)

        return self.__internals__['Green2'][body][(depth, wavenumber)]

    def build_matrices(self, body, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        """Build the matrices of Green coefficients.
        """
        S = np.zeros((self.nb_faces, body.nb_faces), dtype=np.complex64)
        V = np.zeros((self.nb_faces, body.nb_faces), dtype=np.complex64)

        S0, V0 = self._build_matrices_0(body)
        S += S0
        V += V0

        if free_surface < np.infty:

            S1, V1 = self._build_matrices_1(body, free_surface, sea_bottom)
            S += S1
            V += V1

            S2, V2 = self._build_matrices_2(body, free_surface, sea_bottom, wavenumber)
            S += S2
            V += V2

        return S, V

    def get_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Use Meshmagick mesh clipper to remove the part of the mesh above the
        free surface and the part of the mesh below the sea bottom.
        """
        from meshmagick.geometry import Plane
        from meshmagick.mesh_clipper import MeshClipper

        clipped_mesh = MeshClipper(self,
                                   plane=Plane(normal=(0.0, 0.0, 1.0),
                                               scalar=free_surface)).clipped_mesh

        if sea_bottom > -np.infty:
            clipped_mesh = MeshClipper(clipped_mesh,
                                       plane=Plane(normal=(0.0, 0.0, -1.0),
                                                   scalar=-sea_bottom)).clipped_mesh

        clipped_mesh.remove_unused_vertices()

        return FloattingBody(clipped_mesh.vertices, clipped_mesh.faces)
