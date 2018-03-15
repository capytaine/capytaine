#!/usr/bin/env python
# coding: utf-8
"""Special bodies using their symmetries to speed up the computations.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging

import numpy as np

from meshmagick.geometry import Plane

from capytaine.bodies import FloatingBody
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
        half : FloatingBody
            a FloatingBody instance describing half of the body
        plane : Plane
            the symmetry plane across which the half body is mirrored
        name : string, optional
            a name for the body
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


class TranslationalSymmetry(_SymmetricBody):
    """A body composed of a pattern repeated and translated."""

    def __init__(self, body_slice, translation, nb_repetitions=1, name=None):
        """Initialize the body.

        Parameters
        ----------
        body_slice : FloatingBody
            the pattern that will be repeated to form the whole body
        translation : array(3)
            the vector of the translation
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : string, optional
            a name for the body
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


class AxialSymmetry(_SymmetricBody):
    """A body composed of a pattern rotated around a vertical axis."""

    def __init__(self, body_slice, point_on_rotation_axis=np.zeros(3), nb_repetitions=1, name=None):
        """Initialize the body.

        Parameters
        ----------
        body_slice : FloatingBody
            the pattern that will be repeated to form the whole body
        point_on_rotation_axis : array(3)
            one point on the rotation axis. The axis is supposed to be vertical.
        nb_repetitions : int, optional
            the number of repetitions of the pattern (excluding the original one, default: 1)
        name : string, optional
            a name for the body
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
