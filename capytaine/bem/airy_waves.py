"""Computing the potential and velocity of Airy wave."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import numpy as np
from capytaine.tools.lists_of_points import _normalize_points, _normalize_free_surface_points

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
    h = pb.water_depth
    beta = pb.encounter_wave_direction
    wbar = x * np.cos(beta) + y * np.sin(beta)

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
    h = pb.water_depth
    beta = pb.encounter_wave_direction

    wbar = x * np.cos(beta) + y * np.sin(beta)

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
