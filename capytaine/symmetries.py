#!/usr/bin/env python
# coding: utf-8

import numpy as np

from meshmagick.geometry import Plane

from capytaine.bodies import FloattingBody

yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)

class PlanarSymmetry(FloattingBody):

    def __init__(self, half, plane):
        assert isinstance(half, FloattingBody)
        assert isinstance(plane, Plane)

        self.half = half
        self.other_half = self.half.copy()
        self.other_half.mirror(plane)

        self.symmetry_plane = plane

        # self.half.vertices - 2*np.outer(
        #     np.dot(self.half.vertices, self.symmetry_plane.normal) - self.symmetry_plane.c,
        #     self.symmetry_plane.normal
        # )

        self.dofs = {}
        for name, dof in half.dofs:
            self.dofs['Mirrored_' + name] = np.concatenate([dof, -dof])

    @property
    def nb_vertices(self):
        return 2*self.half.nb_vertices

    @property
    def nb_faces(self):
        return 2*self.half.nb_faces

    @property
    def vertices(self):
        return np.concatenate([self.half.vertices, self.other_half.vertices])

    @property
    def faces(self):
        return np.concatenate([self.half.faces, self.other_half.faces + self.half.nb_vertices])

    @property
    def faces_normals(self):
        return np.concatenate([self.half.faces_normals, self.other_half.faces_normals])

    @property
    def faces_areas(self):
        return np.concatenate([self.half.faces_areas, self.other_half.faces_areas])

    @property
    def faces_centers(self):
        return np.concatenate([self.half.faces_centers, self.other_half.faces_centers])

    @property
    def faces_radiuses(self):
        return np.concatenate([self.half.faces_radiuses, self.other_half.faces_radiuses])

    def mirror(self, plane):
        self.half.mirror(plane)
        self.other_half.mirror(plane)
        return

    def build_matrices(self, body, force_full_computation=False, **kwargs):

        if body == self and not force_full_computation:
            Sh, Vh = self.half.build_matrices(self.half, **kwargs)
            Soh, Voh = self.half.build_matrices(self.other_half, **kwargs)

            S = np.concatenate(
                    [np.concatenate([Sh, Soh], axis=1),
                        np.concatenate([Soh, Sh], axis=1)],
                    axis=0)
            V = np.concatenate(
                    [np.concatenate([Vh, Voh], axis=1),
                        np.concatenate([Voh, Vh], axis=1)],
                    axis=0)

        else:
            S1, V1 = self.half.build_matrices(body, **kwargs)
            S2, V2 = self.other_half.build_matrices(body, **kwargs)

            S = np.concatenate(
                    [S1, S2],
                    axis=0)
            V = np.concatenate(
                    [V1, V2],
                    axis=0)

        return S, V


class TranslationSymmetry(FloattingBody):

    def __init__(self, body_slice, translation, nb_repetitions=1):

        assert isinstance(body_slice, FloattingBody)
        assert translation.shape == (3,)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        self.slices = []
        for i in range(nb_repetitions+1):
            new_slice = body_slice.copy()
            new_slice.translate(i*translation)
            self.slices.append(new_slice)

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["Translated_" + name] = np.concatenate([dof]*nb_repetitions)

    @property
    def nb_slices(self):
        return len(self.slices)

    @property
    def nb_vertices(self):
        return self.slices[0].nb_vertices * self.nb_slices

    @property
    def nb_faces(self):
        return self.slices[0].nb_faces * self.nb_slices

    @property
    def vertices(self):
        return np.concatenate([body_slice.vertices for body_slice in self.slices])

    @property
    def faces(self):
        return np.concatenate([
            body_slice.faces + i*self.slices[0].nb_vertices
            for i, body_slice in enumerate(self.slices)
        ])

    @property
    def faces_normals(self):
        return np.concatenate([body_slice.faces_normals for body_slice in self.slices])

    @property
    def faces_areas(self):
        return np.concatenate([body_slice.faces_areas for body_slice in self.slices])

    @property
    def faces_centers(self):
        return np.concatenate([body_slice.faces_centers for body_slice in self.slices])

    @property
    def faces_radiuses(self):
        return np.concatenate([body_slice.faces_radiuses for body_slice in self.slices])

    def mirror(self, plane):
        for body_slice in self.slices:
            body_slice.mirror(plane)
        return

    def build_matrices(self, body, force_full_computation=False, **kwargs):
        """Compute the influence matrix of `self` on `body`.

        force_full_computation (boolean): if True, do not use the symmetry (for debugging).
        """

        if body == self and not force_full_computation:
            Ss, Vs = [], []
            for body_slice in self.slices:
                Si, Vi = self.slices[0].build_matrices(body_slice, **kwargs)
                Ss.append(Si)
                Vs.append(Vi)

            Ss = Ss[:0:-1] + Ss

            S = np.concatenate(
                    [
                        np.concatenate(
                            Ss[i:i+self.nb_slices],
                            axis=1)
                        for i in range(self.nb_slices-1, -1, -1)
                        ],
                    axis=0)

            Vs = Vs[:0:-1] + Vs
            V = np.concatenate(
                    [
                        np.concatenate(
                            Vs[i:i+self.nb_slices],
                            axis=1)
                        for i in range(self.nb_slices-1, -1, -1)
                        ],
                    axis=0)

        else:

            Slist, Vlist = [], []
            for body_slice in self.slices:
                Si, Vi = body_slice.build_matrices(body, **kwargs)
                Slist.append(Si)
                Vlist.append(Vi)

            S = np.concatenate(Slist)
            V = np.concatenate(Vlist)

        return S, V
