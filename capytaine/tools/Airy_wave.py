#!/usr/bin/env python
# coding: utf-8
"""
Compute the potential and velocity of Airy wave.
"""

import numpy as np


def Airy_wave_potential(X, pb):
    """Compute the potential for Airy waves at a given point (or array of points).

    Parameters
    ----------
    X: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (1) or (N x 1)
        The potential
    """
    x, y, z = X.T
    k = pb.wavenumber
    h = pb.depth
    wbar = x*np.cos(pb.angle) + y*np.sin(pb.angle)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        # sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        # sih = np.exp(k*z)

    return -1j*pb.g/pb.omega * cih * np.exp(1j * k * wbar)


def Airy_wave_velocity(X, pb):
    """Compute the fluid velocity for Airy waves at a given point (or array of points).

    Parameters
    ----------
    X: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (3) or (N x 3)
        the velocity vectors
    """
    x, y, z = X.T
    k = pb.wavenumber
    h = pb.depth

    wbar = x*np.cos(pb.angle) + y*np.sin(pb.angle)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        sih = np.exp(k*z)

    v = pb.g*k/pb.omega * \
        np.exp(1j * k * wbar) * \
        np.array([np.cos(pb.angle)*cih, np.sin(pb.angle)*cih, -1j*sih])

    return v.T


def Froude_Krylov_force(problem):
    pressure = -1j * problem.omega * problem.rho * Airy_wave_potential(problem.body.mesh.faces_centers, problem)
    forces = {}
    for dof in problem.influenced_dofs:
        influenced_dof = np.sum(problem.body.dofs[dof] * problem.body.mesh.faces_normals, axis=1)
        forces[dof] = pressure @ (influenced_dof * problem.body.mesh.faces_areas)
    return forces

