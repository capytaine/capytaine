#!/usr/bin/env python
# coding: utf-8

import numpy as np
from capytaine.bodies import FloattingBody
from meshmagick.geometry import Plane

yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)

def reflect_vector(vector, plane):
    return vector - 2*np.dot(vector, plane.normal)*plane.normal

def reflect_point(point, plane):
    return plane.c*plane.normal  + reflect_vector(point - plane.c*plane.normal, plane)

class PlanarSymmetry(FloattingBody):
    """
    """

    def __init__(self, half, plane):
        assert isinstance(half, FloattingBody)
        assert isinstance(plane, Plane)

        self.half = half
        self.symmetry_plane = plane
        self.dofs = {}

    @property
    def nb_vertices(self):
        return 2*self.half.nb_vertices

    @property
    def nb_faces(self):
        return 2*self.half.nb_faces

    @property
    def vertices(self):
        return np.r_[
                self.half.vertices,
                np.asarray(
                    [reflect_point(point, self.symmetry_plane) for point in self.half.vertices]
                    )
                ]

    @property
    def faces(self):
        return np.r_[
                self.half.faces,
                np.fliplr(self.half.faces) + self.half.nb_vertices
                ]

    @property
    def faces_normals(self):
        """Get the array of faces normals of the mesh
        """
        return np.r_[
                self.half.faces_normals,
                np.asarray(
                    [reflect_vector(vector, self.symmetry_plane) for vector in self.half.faces_normals]
                    )
                ]

    @property
    def faces_areas(self):
        return np.r_[
                self.half.faces_areas,
                self.half.faces_areas
                ]

    @property
    def faces_centers(self):
        return np.r_[
                self.half.faces_centers,
                np.asarray(
                    [reflect_point(point, self.symmetry_plane) for point in self.half.faces_centers]
                    )
                ]

    @property
    def faces_radiuses(self):
        return np.r_[
                self.half.faces_radiuses,
                self.half.faces_radiuses,
                ]

    def build_matrices(self, body, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):

        other_half = self.half.copy()
        other_half.mirror(self.symmetry_plane)

        if body == self:
            Sh, Vh = self.half.build_matrices(self.half, free_surface, sea_bottom, wavenumber)
            Soh, Voh = self.half.build_matrices(other_half, free_surface, sea_bottom, wavenumber)
            S = np.r_[
                np.c_[Sh, Soh],
                np.c_[Soh, Sh]
            ]
            V = np.r_[
                np.c_[Vh, Voh],
                np.c_[Voh, Vh]
            ]

        else:
            S1, V1 = self.half.build_matrices(body, free_surface, sea_bottom, wavenumber)
            S2, V2 = other_half.build_matrices(body, free_surface, sea_bottom, wavenumber)

            S = np.r_[S1, S2]
            V = np.r_[V1, V2]

        return S, V
