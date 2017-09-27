#!/usr/bin/env python
# coding: utf-8
"""Special bodies using symmetries to speed up the computations.
"""

import logging

import numpy as np

from meshmagick.geometry import Plane

from capytaine.bodies import FloatingBody


LOG = logging.getLogger(__name__)


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)


class ReflectionSymmetry(FloatingBody):
    """A body composed of two symmetrical halves."""

    def __init__(self, half, plane):
        """Initialize the body.
        
        Parameters
        ----------
        half: FloatingBody
            a FloatingBody instance describing half of the body
        plane: Plane
            the symmetry plane across which the half body is mirrored
        """
        assert isinstance(half, FloatingBody)
        assert isinstance(plane, Plane)

        self.half = half
        self.half.nb_matrices_to_keep *= 2

        self._name = "mirrored_" + half.name
        LOG.info(f"New mirror symmetry: {self.name}.")

        self.other_half = self.half.copy()
        self.other_half.mirror(plane)
        self.other_half.name = "other_half_of_" + self.name

        self.symmetry_plane = plane

        self.dofs = {}
        for name, dof in half.dofs.items():
            self.dofs['mirrored_' + name] = np.concatenate([dof, dof])


    @property
    def nb_matrices_to_keep(self):
        return self.half.nb_matrices_to_keep

    @nb_matrices_to_keep.setter
    def nb_matrices_to_keep(self, value):
        self.half.nb_matrices_to_keep = value
        self.other_half.nb_matrices_to_keep = value

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

    def translate(self, vector):
        self.half.translate(vector)
        self.other_half.translate(vector)
        return

    def build_matrices(self, body, force_full_computation=False, **kwargs):
        """Return the influence matrices of self on body."""
        if body == self and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on itself using mirror symmetry.")

            S = np.empty((self.nb_faces, body.nb_faces), dtype=np.complex64)
            V = np.empty((self.nb_faces, body.nb_faces), dtype=np.complex64)

            # Indices ranges of the four quarters of the matrix
            top_left     = (slice(None, self.nb_faces//2), slice(None, self.nb_faces//2))
            top_right    = (slice(None, self.nb_faces//2), slice(self.nb_faces//2, None))
            bottom_left  = (slice(self.nb_faces//2, None), slice(None, self.nb_faces//2))
            bottom_right = (slice(self.nb_faces//2, None), slice(self.nb_faces//2, None))

            # Evaluation of two of the quarters
            S[top_left], V[top_left] = self.half.build_matrices(self.half, **kwargs)
            S[top_right], V[top_right] = self.half.build_matrices(self.other_half, **kwargs)

            # Copy the values in the two other quarters
            S[bottom_left], V[bottom_left] = S[top_right], V[top_right]
            S[bottom_right], V[bottom_right] = S[top_left], V[top_left]

        else:
            # Do not use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on {body}.")

            S = np.empty((self.nb_faces, body.nb_faces), dtype=np.complex64)
            V = np.empty((self.nb_faces, body.nb_faces), dtype=np.complex64)

            top    = (slice(None, self.nb_faces//2), slice(None, None))
            bottom = (slice(self.nb_faces//2, None), slice(None, None))

            S[top], V[top] = self.half.build_matrices(body, **kwargs)
            S[bottom], V[bottom] = self.other_half.build_matrices(body, **kwargs)

        return S, V


class TranslationalSymmetry(FloatingBody):
    """A body composed of a pattern repeated and translated."""

    def __init__(self, body_slice, translation, nb_repetitions=1):
        """Initialize the body.
        
        Parameters
        ----------
        body_slice: FloatingBody
            the pattern that will be repeated to form the whole body
        translation: array(3)
            the vector of the translation
        nb_repetitions: int
            the number of repetitions of the pattern (excluding the original one)
        """
        assert isinstance(body_slice, FloatingBody)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        translation = np.asarray(translation)
        assert translation.shape == (3,)

        self._name = "translated_" + body_slice.name
        LOG.info(f"New translation symmetry: {self.name}.")

        body_slice.nb_matrices_to_keep *= nb_repetitions
        self.slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy()
            new_slice.translate(i*translation)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            new_slice.name = f"repetition_{i}_of_{body_slice.name}"
            self.slices.append(new_slice)

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["translated_" + name] = np.concatenate([dof]*nb_repetitions)

    @property
    def nb_matrices_to_keep(self):
        return self.slices[0].nb_matrices_to_keep

    @nb_matrices_to_keep.setter
    def nb_matrices_to_keep(self, value):
        for body_slice in self.slices:
            body_slice.nb_matrices_to_keep = value

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

        Parameters
        ----------
        body: FloatingBody
            the body interacting with `self`
        force_full_computation: boolean
            if True, do not use the symmetry (for debugging).
        """

        if body == self and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on itself using translation symmetry.")

            # Compute interactions of all slices with the first one
            Slist, Vlist = [], []
            for body_slice in self.slices:
                Si, Vi = self.slices[0].build_matrices(body_slice, **kwargs)
                Slist.append(Si)
                Vlist.append(Vi)

            # Concatenate elements of the list to build the matrix
            Slist = Slist[:0:-1] + Slist
            S = np.concatenate(
                    [
                        np.concatenate(
                            Slist[i:i+self.nb_slices],
                            axis=1)
                        for i in range(self.nb_slices-1, -1, -1)
                        ],
                    axis=0)

            Vlist = Vlist[:0:-1] + Vlist
            V = np.concatenate(
                    [
                        np.concatenate(
                            Vlist[i:i+self.nb_slices],
                            axis=1)
                        for i in range(self.nb_slices-1, -1, -1)
                        ],
                    axis=0)

        else:
            # Do not use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on {body}.")

            Slist, Vlist = [], []
            for body_slice in self.slices:
                Si, Vi = body_slice.build_matrices(body, **kwargs)
                Slist.append(Si)
                Vlist.append(Vi)

            S = np.concatenate(Slist)
            V = np.concatenate(Vlist)

        return S, V
