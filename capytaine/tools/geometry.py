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
    def __init__(self, normal=(0., 0., 1.), scalar=0., name=None):

        normal = np.asarray(normal, dtype=np.float)

        self._normal = normal / np.linalg.norm(normal)
        self._scalar = float(scalar)

        self.name = str(name)

    def __str__(self):
        str_repr = "Plane{normal=[%f, %f, %f], scalar=%f}" % \
                   (self._normal[0], self._normal[1], self._normal[2], self._scalar)
        return str_repr

    def __eq__(self, other):
        if isinstance(other, Plane):
            return (np.isclose(self.c, other.c, atol=1e-5) and
                    np.isclose(self.normal @ other.normal, 1.0, atol=1e-5))
        else:
            raise NotImplemented()

    @property
    def normal(self):
        """Get the plane's normal"""
        return self._normal

    @normal.setter
    def normal(self, value):
        """Set the plane's normal"""
        value = np.asarray(value, dtype=np.float)
        self._normal = value / np.linalg.norm(value)

    @property
    def c(self):
        """Get the plane's scalar parameter"""
        return self._scalar

    @c.setter
    def c(self, value):
        """Set the scalar parameter of the plane equation"""
        self._scalar = float(value)

    def rotate(self, axis, angle):
        rot_matrix = axis.rotation_matrix(angle)
        self.normal = rot_matrix @ self.normal

    def set_normal_from_angles(self, theta_x, theta_y):
        """Set the normal orientation given angles theta_x and theta_y.

        Parameters
        ----------
        theta_x : float
            Angle around Ox (rad)
        theta_y : float
            Angle around Oy (rad)
        """

        theta = np.sqrt(theta_x * theta_x + theta_y * theta_y)

        if theta == 0.:
            self.normal[:] = [0., 0., 1.]
        else:
            stheta_theta = np.sin(theta) / theta
            ctheta = np.cos(theta)
            self.normal[:] = np.array([stheta_theta * theta_y,
                                       -stheta_theta * theta_x,
                                       ctheta])

    def get_normal_orientation_wrt_z(self):
        """Returns the angles theta_x and theta_y giving the orientation of the plane normal"""

        nx, ny, nz = self.normal
        stheta = np.sqrt(nx*nx + ny*ny)
        ctheta = nz
        theta_x = theta_y = 0.
        if stheta == 0.:
            if nz == 1.:
                theta_x = theta_y = 0.
            elif nz == -1.:
                theta_x = np.pi
                theta_y = 0.
        else:
            theta = atan2(stheta, ctheta)
            theta_stheta = theta / stheta

            theta_x = -theta_stheta * ny
            theta_y = theta_stheta * nx

        return theta_x, theta_y

    def set_plane_parameters(self, scalar, theta_x, theta_y):
        """
        Updates the plane parameters (normal and scalar parameter) given scalar and angles.

        Parameters
        ----------
        scalar : float
            Plane scalar parameter (m)
        theta_x : float
            Normal angle around Ox (rad)
        theta_y : float
            Normal angle around Oy (rad)
        """

        self.rotate_normal(theta_x, theta_y)
        ctheta = np.cos(np.sqrt(theta_x * theta_x + theta_y * theta_y))
        self._scalar = self._scalar * ctheta + scalar

    def get_point_dist_wrt_plane(self, points):
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

        return np.dot(points, self._normal) - self._scalar

    def flip_normal(self):
        """
        Flips the Normal of the plane
        """
        self.normal *= -1
        theta_x, theta_y = self.get_normal_orientation_wrt_z()

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
        t = (p0n - self._scalar) / (p0n - p1n)
        if t < 0. or t > 1.:
            raise RuntimeError('Intersection is outside the edge')
        return (1-t) * p0 + t * p1

    def orthogonal_projection_on_plane(self, points):
        """
        Returns the coordinates of the orthogonal projection of points

        Parameters
        ----------
        points : ndarray
            Coordinates of the points to be projected

        Returns
        -------
        projected_points : ndarray
            Coordinates of the projection points
        """
        # TODO: passer en vectoriel
        projected_points = np.zeros_like(points)
        for point, projected_point in zip(points, projected_points):
            projected_point[:] = point - self.get_point_dist_wrt_plane(point) * self.normal

        return projected_points

    def get_origin(self):
        """Get the coordinates of the plane's origin"""
        return self.c * self.normal


# Useful aliases
yOz_Plane = Plane(normal=(1.0, 0.0, 0.0), scalar=0.0)
xOz_Plane = Plane(normal=(0.0, 1.0, 0.0), scalar=0.0)
xOy_Plane = Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)

def _cardan(phi, theta):
    """
    Computes the rotation matrix corresponding to angles phi (roll) and theta (pitch) using Cardan angles convention

    Parameters
    ----------
    phi : float
        Roll angle (rad)
    theta : float
        Pitch angel (rad)

    Returns
    -------
    R : ndarray
        Rotation matrix
    """

    # Rotation matrix
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    rot_e0 = np.zeros((3, 3), dtype=float)
    rot_e0[0] = [ctheta, 0., -stheta]
    rot_e0[1] = [sphi*stheta, cphi, sphi*ctheta]
    rot_e0[2] = [cphi*stheta, -sphi, cphi*ctheta]

    return rot_e0
