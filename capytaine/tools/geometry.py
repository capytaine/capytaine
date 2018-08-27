#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Freely adapted from meshmagick.
"""

from abc import ABC, abstractmethod
from math import atan2

import numpy as np


def inplace_transformation(inplace_function):
    """Decorator for methods transforming 3D objects:
    * Add the optional argument `inplace` to return a new object instead of doing the transformation in place.
    * If the object has properties cached in an "__internals__" dict, they are deleted.
    """
    def enhanced_inplace_function(self, *args, inplace=True, name=None, **kwargs):
        if not inplace:
            object3d = self.copy(name=name)
        else:
            object3d = self
        inplace_function(object3d, *args, **kwargs)
        if hasattr(object3d, '__internals__'):
            object3d.__internals__.clear()
        return object3d
    return enhanced_inplace_function


class Abstract3DObject(ABC):
    """Abstract class for 3d objects that can be transformed in 3d.
    The child class have to define mirror, rotate and translate,
    then more routines such as translate_x are automatically defined."""

    @abstractmethod
    def translate(self, vector):
        pass

    @abstractmethod
    def rotate(self, axis, angle):
        pass

    @abstractmethod
    def mirror(self, plane):
        pass

    @inplace_transformation
    def translate_x(self, tx):
        return self.translate((tx, 0., 0.))

    @inplace_transformation
    def translate_y(self, ty):
        return self.translate((0., ty, 0.))

    @inplace_transformation
    def translate_z(self, tz):
        return self.translate((0., 0., tz))

    @inplace_transformation
    def rotate_x(self, thetax):
        return self.rotate(Ox_axis, thetax)

    @inplace_transformation
    def rotate_y(self, thetay):
        return self.rotate(Oy_axis, thetay)

    @inplace_transformation
    def rotate_z(self, thetaz):
        return self.rotate(Oz_axis, thetaz)

    @inplace_transformation
    def rotate_angles(self, angles):
        thetax, thetay, thetaz = angles
        self.rotate(Ox_axis, thetax)
        self.rotate(Oy_axis, thetay)
        self.rotate(Oz_axis, thetaz)
        return self


class Axis(Abstract3DObject):
    def __init__(self, vector, point=(0, 0, 0)):

        assert len(vector) == 3, "Vector of an axis should be given as a 3-ple of values."
        assert len(point) == 3, "Point of an axis should be given as a 3-ple of values."

        self.vector = np.array(vector, np.float)
        self.point = np.array(point, np.float)

    def copy(self, name=None):
        return Axis(vector=self.vector, point=self.point)

    def __contains__(self, point):
        assert len(point) == 3, "Points should be given as a 3-ple of values."
        point = np.asarray(point, dtype=np.float)
        return np.isclose(np.linalg.norm(np.cross(point - self.point, self.vector)), 0)

    def is_orthogonal_to(self, item):
        if isinstance(item, Axis):
            raise NotImplementedError
        elif isinstance(item, Plane):
            return np.isclose(np.linalg.norm(np.cross(self.vector, item.normal)), 0.0)
        else:  # The item is supposed to be a vector given as a 3-ple
            return np.isclose(np.linalg.norm(self.vector @ item), 0.0)

    def rotation_matrix(self, theta):
        """Rotation matrix around the vector according to Rodrigues' formula."""
        ux, uy, uz = self.vector
        W = np.array([[0, -uz, uy],
                      [uz, 0, -ux],
                      [-uy, ux, 0]])
        return np.identity(3) + np.sin(theta)*W + 2*np.sin(theta/2)**2 * (W @ W)

    @inplace_transformation
    def translate(self, vector):
        self.point += vector
        return

    @inplace_transformation
    def rotate(self, axis, angle):
        raise NotImplemented

    @inplace_transformation
    def mirror(self, plane):
        raise NotImplemented


Ox_axis = Axis(vector=(1, 0, 0), point=(0, 0, 0))
Oy_axis = Axis(vector=(0, 1, 0), point=(0, 0, 0))
Oz_axis = Axis(vector=(0, 0, 1), point=(0, 0, 0))


class Plane:
    """Class to handle plane geometry.

    A plane is represented by the equation :math:`\\vec{n}.\\vec{x} = c` where :math:`\\vec{n}` is the plane's normal,
    :math:`\\vec{x}` a point in the space and :math:`c` a scalar parameter being the signed distance between the
    reference frame origin and the its otrhogonal projection on the plane.

    Parameters
    ----------
    normal : array_like
        3 component vector of the plane normal
    scalar : float
        The scalar parameter of the plane
    """
    def __init__(self, normal=(0.0, 0.0, 1.0), scalar=0.0, name=None):
        self.normal = normal
        self.c = float(scalar)
        self.name = str(name)

    @property
    def normal(self):
        """Get the plane's normal"""
        return self._normal

    @normal.setter
    def normal(self, value):
        """Set the plane's normal"""
        value = np.asarray(value, dtype=np.float)
        self._normal = value / np.linalg.norm(value)

    def __repr__(self):
        return f"Plane(normal=({self.normal[0]}, {self.normal[1]}, {self.normal[2]}), scalar={self.c})"

    def __eq__(self, other):
        if isinstance(other, Plane):
            return (np.isclose(self.c, other.c, atol=1e-5) and
                    np.isclose(self.normal @ other.normal, 1.0, atol=1e-5))
        else:
            return NotImplemented

    def copy(self, name=None):
        if name is None:
            name = self.name
        return Plane(normal=self.normal, scalar=self.c, name=name)

    @property
    def origin(self):
        """Get the coordinates of the plane's origin"""
        return self.c * self.normal

    def __contains__(self, item):
        if isinstance(item, Axis):
            return item.point in self and item.vector
        assert len(item) == 3, "Points should be given as a 3-ple of values."
        point = np.asarray(point, dtype=np.float)
        return np.isclose(np.linalg.norm(np.cross(point - self.point, self.vector)), 0)

    def is_orthogonal_to(self, other):
        if isinstance(other, Axis):
            return np.isclose(np.linalg.norm(np.cross(self.normal, other.vector)), 0.0)
        elif isinstance(other, Plane):
            return np.isclose(np.linalg.norm(self.normal @ other.normal), 0.0)
        else:  # The other is supposed to be a vector given as a 3-ple
            return np.isclose(np.linalg.norm(np.cross(self.normal, other)), 0.0)

    @inplace_transformation
    def rotate(self, axis, angle):
        if self.c != 0.0 or (0, 0, 0) not in axis:
            raise NotImplementedError
        rot_matrix = axis.rotation_matrix(angle)
        self.normal = rot_matrix @ self.normal
        return self

    @inplace_transformation
    def translate(self, vector):
        self.c = self.c + self.normal @ np.asarray(vector)
        return self

    @inplace_transformation
    def mirror(self, plane):
        raise NotImplementedError

    def distance_to_point(self, points):
        """
        Return the orthogonal distance of points with respect to the plane

        Parameters
        ----------
        points : ndarray
            Array of points coordinates

        Returns
        -------
        dist : ndarray
            Array of distances of points with respect to the plane
        """

        return np.dot(points, self._normal) - self.c

    def get_edge_intersection(self, p0, p1):
        """
        Returns the coordinates of the intersection point between the plane and the edge P0P1.

        Parameters
        ----------
        p0 : ndarray
            Coordinates of point p0
        p1 : ndarray
            Coordinates of point P1

        Returns
        -------
        I : ndarray
            Coordinates of intersection point
        """
        assert len(p0) == 3 and len(p1) == 3

        p0n = np.dot(p0, self.normal)
        p1n = np.dot(p1, self.normal)
        t = (p0n - self.c) / (p0n - p1n)
        if t < 0. or t > 1.:
            raise RuntimeError('Intersection is outside the edge')
        return (1-t) * p0 + t * p1


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)

