#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
From meshmagick.
"""

import numpy as np
from math import atan2

# TODO: voir si on ne peut pas mettre ces fonctions dans un module dedie --> module rotation !!!


def _rotation_matrix(angles):
    angles = np.asarray(angles, dtype=np.float)
    theta = np.linalg.norm(angles)
    if theta == 0.:
        return np.eye(3)

    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    nx, ny, nz = angles/theta
    nxny = nx*ny
    nxnz = nx*nz
    nynz = ny*nz
    nx2 = nx*nx
    ny2 = ny*ny
    nz2 = nz*nz

    return (ctheta*np.eye(3)
           + (1-ctheta) * np.array([[nx2, nxny, nxnz],
                                    [nxny, ny2, nynz],
                                    [nxnz, nynz, nz2]])
           + stheta * np.array([[0., -nz, ny],
                                [nz, 0., -nx],
                                [-ny, nx, 0.]])
            )


def _rodrigues(thetax, thetay):
    """
    Computes the rotation matrix corresponding to angles thetax and thetay using the Olinde-Rodrigues formula

    Parameters
    ----------
    thetax : float
        Angle around Ox axe (rad)
    thetay
        Angle around Oy axe (rad)
    Returns
    -------
    rot : ndarray
        Rotation matrix
    """

    theta = np.sqrt(thetax*thetax + thetay*thetay)
    if theta == 0.:
        nx = ny = 0.
        ctheta = 1.
        stheta = 0.
    else:
        nx, ny = thetax/theta, thetay/theta
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
    nxny = nx*ny

    # Olinde Rodrigues formulae
    # FIXME: S'assurer qu'on a effectivement pas de Ctheta devant le I3 !! et repercuter sur l'hydrostatique
    rot = ctheta*np.eye(3) \
        + (1-ctheta) * np.array([[nx*nx, nxny, 0.],
                                 [nxny, ny*ny, 0.],
                                 [0., 0., 0.]]) \
        + stheta * np.array([[0., 0.,  ny],
                             [0., 0., -nx],
                             [-ny, nx, 0.]])
    return rot


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


def _get_rotation_matrix(theta_x, theta_y, atype='fixed'):
    """
    Computes rotation matrix using different angle conventions

    Parameters
    ----------
    theta_x : float
        Angle around x (rad)
    theta_y : float
        Angle around y (rad)
    atype : {'fixed', 'cardan'}, optional
        Angle convention to use. Default to 'fixed' (fixed axes)

    Returns
    -------
    ndarray
        Rotation matrix

    """
    if atype == 'fixed':
        rot_matrix = _rodrigues(theta_x, theta_y)
    elif atype == 'cardan':
        rot_matrix = _cardan(theta_x, theta_y)
    else:
        raise AttributeError('Unknown angle convention: %s' % atype)

    return rot_matrix


def _get_axis_angle_from_rotation_matrix(rot_matrix):
    """Returns the angle and unit rotation axis from a rotation matrix"""
    warn('Fonction _get_axis_angle_from_rotation_matrix a verifier !!!')
    theta = np.acos((np.trace(rot_matrix) - 1.) * 0.5)
    direction = (1./(2.*np.sin(theta))) * np.array([rot_matrix[2, 1] - rot_matrix[1, 2],
                                                      rot_matrix[0, 2] - rot_matrix[2, 0],
                                                      rot_matrix[1, 0] - rot_matrix[0, 1]])
    return theta, direction


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

        # Storing rotation matrix (redundant !) to speedup computations
        # Shall be _update in methods !!! --> using decorator ?
        theta_x, theta_y = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(theta_x, theta_y)

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

    def rotate(self, angles):
        rot_matrix = _rotation_matrix(angles)
        self.normal = rot_matrix @ self.normal
        self._rot = rot_matrix @ self._rot

    def rotate_normal(self, theta_x, theta_y):
        """
        Rotates the current plane normal by fixed angles theta_x and theta_y.

        Parameters
        ----------
        theta_x : float
            Angle of rotation around Ox (rad)
        theta_y : float
            Angle of rotation around Oy (rad)
        """
        rot_matrix = _get_rotation_matrix(theta_x, theta_y)
        self.normal = np.dot(rot_matrix, self.normal)

        # updating self._rot
        self._rot = np.dot(rot_matrix, self._rot)

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

            # updating self._rot
            self._rot = np.eye(3)
        else:
            stheta_theta = np.sin(theta) / theta
            ctheta = np.cos(theta)
            self.normal[:] = np.array([stheta_theta * theta_y,
                                       -stheta_theta * theta_x,
                                       ctheta])
            # Updating self._rot
            self._rot = _get_rotation_matrix(theta_x, theta_y)

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
        self._rot = _get_rotation_matrix(theta_x, theta_y)

    def coord_in_plane(self, points):
        """
        Return the coordinates of points in the frame of the plane

        Parameters
        ----------
        points : ndarray
            Array of points coordinates

        Returns
        -------
        output : ndarray
            Array of points coordinates in the frame of the plane
        """

        # TODO: verifier effectivement que si on prend des points se trouvant dans le plan, leurs coordonnees dans le
        #  plan n'ont pas de composante z
        return -self._scalar * self.normal + np.transpose(np.dot(self._rot, points.T))

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

