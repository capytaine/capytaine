#!/usr/bin/env python
# coding: utf-8
"""Computing the potential and velocity of Airy wave."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import numpy as np
from capytaine.meshes import Mesh, CollectionOfMeshes


def _normalize_points(points, keep_mesh=False):
    if isinstance(points, (Mesh, CollectionOfMeshes)):
        if keep_mesh:
            return points, (points.nb_faces,)
        else:
            return points.faces_centers, (points.nb_faces,)

    points = np.asarray(points)

    if points.ndim == 1:  # A single point has been provided
        output_shape = (1,)
        points = points.reshape((1, points.shape[0]))

    elif points.ndim == 2:
        output_shape = (points.shape[0],)

    elif points.ndim > 2:
        # `points` is expected to be the resuls of a meshgrid. Points has shape (d, nx, ny, ...)
        output_shape = points.shape[1:]
        points = points.reshape(points.shape[0], -1).transpose()
        # points is now a (nx*ny*... , d) array

    else:
        raise ValueError("This should not happen.")

    return points, output_shape

def _normalize_free_surface_points(points, keep_mesh=False):
    if keep_mesh and isinstance(points, (Mesh, CollectionOfMeshes)):
        return points, (points.nb_faces,)

    points, output_shape = _normalize_points(points, keep_mesh)

    if points.ndim == 2 and points.shape[1] == 2:  # Only x and y have been provided
        points = np.concatenate([points, np.zeros((points.shape[0], 1))], axis=1)

    return points, output_shape



def airy_waves_potential(points, pb):
    """Compute the potential for Airy waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (1) or (N x 1)
        The potential
    """
    points, output_shape = _normalize_points(points)

    x, y, z = points.T
    k = pb.wavenumber
    h = pb.depth
    wbar = x * np.cos(pb.wave_direction) + y * np.sin(pb.wave_direction)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        # sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        # sih = np.exp(k*z)

    phi = -1j*pb.g/pb.omega * cih * np.exp(1j * k * wbar)
    return phi.reshape(output_shape)


def airy_waves_velocity(points, pb):
    """Compute the fluid velocity for Airy waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (3) or (N x 3)
        the velocity vectors
    """

    points, output_shape = _normalize_points(points)

    x, y, z = points.T
    k = pb.wavenumber
    h = pb.depth

    wbar = x * np.cos(pb.wave_direction) + y * np.sin(pb.wave_direction)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        sih = np.exp(k*z)

    v = pb.g*k/pb.omega * \
        np.exp(1j * k * wbar) * \
        np.array([np.cos(pb.wave_direction) * cih, np.sin(pb.wave_direction) * cih, -1j * sih])

    return v.T.reshape((*output_shape, 3))


def airy_waves_pressure(points, pb):
    return 1j * pb.omega * pb.rho * airy_waves_potential(points, pb)


def froude_krylov_force(pb):
    return pb.body.integrate_pressure(airy_waves_pressure(pb.body.mesh.faces_centers, pb))


def airy_waves_free_surface_elevation(points, pb):
    """Compute the free surface elevation at points of the undisturbed Airy waves

    Parameters
    ----------
    points: array of shape (3) or (N × 3) or (2) or (N × 2)
        coordinates of the points in which to evaluate the potential.
        If only two coordinates are passed, the last one is filled with zeros.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    complex-valued array of shape (1,) or (N,)
        the free surface elevations
    """
    points, output_shape = _normalize_free_surface_points(points)
    return 1j * pb.omega / pb.g * airy_waves_potential(points, pb).reshape(output_shape)


