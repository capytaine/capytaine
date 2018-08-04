#!/usr/bin/env python
# coding: utf-8
"""Compute the potential and velocity of Airy wave."""
# This file is part of "capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import numpy as np


def Airy_wave_potential(X, pb, convention="Nemoh"):
    """Compute the potential for Airy waves at a given point (or array of points).

    Parameters
    ----------
    X: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest
    convention: str, optional
        convention for the incoming wave field. Accepted values: "Nemoh", "WAMIT".

    Returns
    -------
    array of shape (1) or (N x 1)
        The potential
    """
    assert convention.lower() in ["nemoh", "wamit"], \
        "Convention for wave field should be either Nemoh or WAMIT."

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

    if convention.lower() == "wamit":
        return  1j*pb.g/pb.omega * cih * np.exp(-1j * k * wbar)
    else:
        return -1j*pb.g/pb.omega * cih * np.exp(1j * k * wbar)


def Airy_wave_velocity(X, pb, convention="Nemoh"):
    """Compute the fluid velocity for Airy waves at a given point (or array of points).

    Parameters
    ----------
    X: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest
    convention: str, optional
        convention for the incoming wave field. Accepted values: "Nemoh", "WAMIT".

    Returns
    -------
    array of shape (3) or (N x 3)
        the velocity vectors
    """
    assert convention.lower() in ["nemoh", "wamit"], \
        "Convention for wave field should be either Nemoh or WAMIT."

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

    if convention.lower() == "wamit":
        return np.conjugate(v.T)
    else:
        return v.T


def Froude_Krylov_force(problem, convention="Nemoh"):
    pressure = -1j * problem.omega * problem.rho * Airy_wave_potential(problem.body.mesh.faces_centers, problem, convention=convention)
    forces = {}
    for dof in problem.influenced_dofs:
        influenced_dof = np.sum(problem.body.dofs[dof] * problem.body.mesh.faces_normals, axis=1)
        forces[dof] = pressure @ (influenced_dof * problem.body.mesh.faces_areas)
    return forces

