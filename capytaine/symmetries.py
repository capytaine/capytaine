#!/usr/bin/env python
# coding: utf-8
"""Special bodies using their symmetries to speed up the computations."""

import logging

import numpy as np

from meshmagick.geometry import Plane

from capytaine.bodies import FloatingBody
from capytaine.Toeplitz_matrices import BlockToeplitzMatrix, BlockCirculantMatrix
from capytaine.bodies_collection import CollectionOfFloatingBodies


LOG = logging.getLogger(__name__)


class _SymmetricBody(CollectionOfFloatingBodies):
    def __add__(self, body_to_add):
        return CollectionOfFloatingBodies([self, body_to_add])


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)


class ReflectionSymmetry(_SymmetricBody):
    """A body composed of two symmetrical halves."""

    def __init__(self, half, plane, name=None):
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

        self.plane = plane

        other_half = half.copy()
        other_half.mirror(plane)
        other_half.name = "mirror_of_" + half.name

        CollectionOfFloatingBodies.__init__(self, [half, other_half])

        if name is None:
            self.name = f"ReflectionSymmetry({half.name})"
        else:
            self.name = name
        LOG.info(f"New mirror symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in half.dofs.items():
            self.dofs['mirrored_' + name] = np.concatenate([dof, dof])

    def build_matrices(self, other_body, force_full_computation=False, **kwargs):
        """Return the influence matrices of self on other_body."""
        if isinstance(other_body, ReflectionSymmetry) and other_body.plane == self.plane and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            if other_body == self:
                LOG.debug(f"Evaluating matrix of {self.name} on itself using mirror symmetry.")
            else:
                LOG.debug(f"Evaluating matrix of {self.name} on {other_body.name} itself using mirror symmetry.")

            S_a, V_a = self.subbodies[0].build_matrices(other_body.subbodies[0], **kwargs)
            S_b, V_b = self.subbodies[0].build_matrices(other_body.subbodies[1], **kwargs)

            return BlockToeplitzMatrix([S_a, S_b]), BlockToeplitzMatrix([V_a, V_b])

        else:
            return CollectionOfFloatingBodies.build_matrices(self, other_body, **kwargs)


class TranslationalSymmetry(_SymmetricBody):
    """A body composed of a pattern repeated and translated."""

    def __init__(self, body_slice, translation, nb_repetitions=1, name=None):
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

        self.translation = translation

        body_slice.nb_matrices_to_keep *= nb_repetitions+1
        slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy(name=f"repetition_{i}_of_{body_slice.name}")
            new_slice.translate(i*translation)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            slices.append(new_slice)

        CollectionOfFloatingBodies.__init__(self, slices)

        if name is None:
            self.name = f"TranslationSymmetry({body_slice.name})"
        else:
            self.name = name
        LOG.info(f"New translation symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["translated_" + name] = np.concatenate([dof]*nb_repetitions)

    def build_matrices(self, other_body, force_full_computation=False, **kwargs):
        """Compute the influence matrix of `self` on `other_body`.

        Parameters
        ----------
        body: FloatingBody
            the body interacting with `self`
        force_full_computation: boolean
            if True, do not use the symmetry (for debugging).
        """

        if (isinstance(other_body, TranslationalSymmetry)
                and np.allclose(other_body.translation, self.translation)
                and other_body.nb_subbodies == self.nb_subbodies
                and not force_full_computation):
            # Use symmetry to speed up the evaluation of the matrix
            if other_body == self:
                LOG.debug(f"Evaluating matrix of {self.name} on itself using translation symmetry.")
            else:
                LOG.debug(f"Evaluating matrix of {self.name} on {other_body.name} itself using translation symmetry.")

            S_list, V_list = [], []
            for body in other_body.subbodies:
                S, V = self.subbodies[0].build_matrices(body, **kwargs)
                S_list.append(S)
                V_list.append(V)
            return BlockToeplitzMatrix(S_list), BlockToeplitzMatrix(V_list)

        else:
            return CollectionOfFloatingBodies.build_matrices(self, other_body, **kwargs)


class AxialSymmetry(_SymmetricBody):
    """A body composed of a pattern rotated around a vertical axis."""

    def __init__(self, body_slice, point_on_rotation_axis=np.zeros(3), nb_repetitions=1, name=None):
        """Initialize the body.

        Parameters
        ----------
        body_slice: FloatingBody
            the pattern that will be repeated to form the whole body
        point_on_rotation_axis: array(3)
            one point on the rotation axis. The axis is supposed to be vertical.
        nb_repetitions: int
            the number of repetitions of the pattern (excluding the original one)
        """
        assert isinstance(body_slice, FloatingBody)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        point_on_rotation_axis = np.asarray(point_on_rotation_axis)
        assert point_on_rotation_axis.shape == (3,)

        body_slice.nb_matrices_to_keep *= nb_repetitions+1
        slices = [body_slice]
        for i in range(1, nb_repetitions+1):
            new_slice = body_slice.copy(name=f"rotation_{i}_of_{body_slice.name}")
            new_slice.translate(-point_on_rotation_axis)
            new_slice.rotate_z(2*i*np.pi/(nb_repetitions+1))
            new_slice.translate(point_on_rotation_axis)
            new_slice.nb_matrices_to_keep *= nb_repetitions+1
            slices.append(new_slice)

        CollectionOfFloatingBodies.__init__(self, slices)

        if name is None:
            self.name = f"AxialSymmetry({body_slice.name})"
        else:
            self.name = name
        LOG.info(f"New rotation symmetric body: {self.name}.")

        self.dofs = {}
        for name, dof in body_slice.dofs.items():
            self.dofs["rotated_" + name] = np.concatenate([dof]*nb_repetitions)

    def build_matrices(self, other_body, force_full_computation=False, **kwargs):
        """Compute the influence matrix of `self` on `other_body`.

        Parameters
        ----------
        other_body: FloatingBody
            the body interacting with `self`
        force_full_computation: boolean
            if True, do not use the symmetry (for debugging).
        """

        if other_body == self and not force_full_computation:
            # Use symmetry to speed up the evaluation of the matrix
            LOG.debug(f"Evaluating matrix of {self.name} on itself using rotation symmetry.")

            S_list, V_list = [], []
            for body in self.subbodies[:self.nb_subbodies//2+1]:
                S, V = self.subbodies[0].build_matrices(body, **kwargs)
                S_list.append(S)
                V_list.append(V)
            if self.nb_subbodies % 2 == 0:
                return BlockToeplitzMatrix(S_list + S_list[-2:0:-1]), BlockToeplitzMatrix(V_list + V_list[-2:0:-1])
            else:
                return BlockToeplitzMatrix(S_list + S_list[-1:0:-1]), BlockToeplitzMatrix(V_list + V_list[-1:0:-1])

        else:
            return CollectionOfFloatingBodies.build_matrices(self, other_body, **kwargs)
