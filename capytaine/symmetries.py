#!/usr/bin/env python
# coding: utf-8
"""Special bodies using symmetries to speed up the computations."""

import logging
from itertools import chain, accumulate

import numpy as np

from meshmagick.geometry import Plane

from capytaine.bodies import FloatingBody


LOG = logging.getLogger(__name__)


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)


class CollectionOfFloatingBodies(FloatingBody):
    """A body composed of several floating bodies."""

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
        nb_faces = accumulate(chain([0], (body.nb_vertices for body in self.subbodies[:-1])))
        total_nb_faces = sum(nb_faces)
        for nbf, body in zip(nb_faces, bodies):
            for name, dof in body.dofs.items():
                self.dofs['_'.join([body.name, name])] = np.r_[
                    np.zeros(nbf),
                    dof,
                    np.zeros(total_nb_faces - len(dof) - nbf),
                ]

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

    def translate(self, vector):
        for body in self.subbodies:
            body.translate(vector)
        return

    def build_matrices(self, other_body, **kwargs):
        LOG.debug(f"Evaluating matrix of {self.name} on {other_body.name}.")

        S = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
        V = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

        nb_faces = list(accumulate(chain([0], (body.nb_faces for body in self.subbodies))))
        for (i, j), body in zip(zip(nb_faces, nb_faces[1:]), self.subbodies):
            matrix_slice = (slice(i, j), slice(None, None))
            S[matrix_slice], V[matrix_slice] = body.build_matrices(other_body, **kwargs)

        return S, V


class ReflectionSymmetry(CollectionOfFloatingBodies):
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

        half.nb_matrices_to_keep *= 2

        other_half = half.copy()
        other_half.mirror(plane)
        other_half.name = "mirror_of_" + half.name

        CollectionOfFloatingBodies.__init__(self, [half, other_half])

        self.name = "mirrored_" + half.name
        LOG.info(f"New mirror symmetry: {self.name}.")

        self.dofs = {}
        for name, dof in half.dofs.items():
            self.dofs['mirrored_' + name] = np.concatenate([dof, dof])

    def build_matrices(self, other_body, force_full_computation=False, **kwargs):
        """Return the influence matrices of self on other_body."""
        if other_body == self and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on itself using mirror symmetry.")

            S = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
            V = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

            # Indices ranges of the four quarters of the matrix
            top_left     = (slice(None, self.nb_faces//2), slice(None, self.nb_faces//2))
            top_right    = (slice(None, self.nb_faces//2), slice(self.nb_faces//2, None))
            bottom_left  = (slice(self.nb_faces//2, None), slice(None, self.nb_faces//2))
            bottom_right = (slice(self.nb_faces//2, None), slice(self.nb_faces//2, None))

            # Evaluation of two of the quarters
            S[top_left], V[top_left] = self.subbodies[0].build_matrices(self.subbodies[0], **kwargs)
            S[top_right], V[top_right] = self.subbodies[0].build_matrices(self.subbodies[1], **kwargs)

            # Copy the values in the two other quarters
            S[bottom_left], V[bottom_left] = S[top_right], V[top_right]
            S[bottom_right], V[bottom_right] = S[top_left], V[top_left]

            return S, V

        else:
            return CollectionOfFloatingBodies.build_matrices(self, other_body, **kwargs)


class TranslationalSymmetry(CollectionOfFloatingBodies):
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

        body_slice.nb_matrices_to_keep *= nb_repetitions
        slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy()
            new_slice.translate(i*translation)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            new_slice.name = f"repetition_{i}_of_{body_slice.name}"
            slices.append(new_slice)

        CollectionOfFloatingBodies.__init__(self, slices)

        self.name = "repeated_" + body_slice.name
        LOG.info(f"New translation symmetry: {self.name}.")

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["translated_" + name] = np.concatenate([dof]*nb_repetitions)

    def build_matrices(self, other_body, force_full_computation=False, **kwargs):
        """Compute the influence matrix of `self` on `body`.

        Parameters
        ----------
        body: FloatingBody
            the body interacting with `self`
        force_full_computation: boolean
            if True, do not use the symmetry (for debugging).
        """

        if other_body == self and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on itself using translation symmetry.")

            S = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)
            V = np.empty((self.nb_faces, other_body.nb_faces), dtype=np.complex64)

            # Compute the indices of each block composing the block matrix
            nb_faces = list(accumulate(chain([0], (body.nb_faces for body in self.subbodies))))
            matrix_blocks = []
            for i, j in zip(nb_faces, nb_faces[1:]):
                line_of_blocks = []
                for k, l in zip(nb_faces, nb_faces[1:]):
                    line_of_blocks.append((slice(i, j), slice(k, l)))
                matrix_blocks.append(line_of_blocks)

            # Compute the values of the first column of blocks
            for matrix_block, body in zip(matrix_blocks[0], self.subbodies):
                S[matrix_block], V[matrix_block] = self.subbodies[0].build_matrices(body, **kwargs)

            # Copy in the rest of the matrix
            for i in range(1, len(matrix_blocks)):
                for j in range(len(matrix_blocks[i])):
                    S[matrix_blocks[i][j]] = S[matrix_blocks[0][abs(j-i)]]
                    V[matrix_blocks[i][j]] = V[matrix_blocks[0][abs(j-i)]]

            return S, V

        else:
            return CollectionOfFloatingBodies.build_matrices(self, body, **kwargs)


