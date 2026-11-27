# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Computing the potential and velocity of Airy wave."""

import numpy as np
from capytaine.tools.lists_of_points import _normalize_points, _normalize_free_surface_points

def airy_waves_potential(points, problem):
    """Compute the potential for Airy waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    problem: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (1) or (N x 1)
        The potential
    """
    points, output_shape = _normalize_points(points)

    if float(problem.wavenumber) == 0.0 or float(problem.wavenumber) == np.inf:
        return np.nan * np.ones(output_shape)

    x, y, z = points.T
    k = problem.wavenumber
    h = problem.water_depth
    beta = problem.encounter_wave_direction
    wbar = x * np.cos(beta) + y * np.sin(beta)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        # sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        # sih = np.exp(k*z)

    phi = -1j*problem.g/problem.omega * cih * np.exp(1j * k * wbar)
    return phi.reshape(output_shape)


def airy_waves_velocity(points, problem):
    """Compute the fluid velocity for Airy waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    problem: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (3) or (N x 3)
        the velocity vectors
    """
    points, output_shape = _normalize_points(points)

    if float(problem.wavenumber) == 0.0 or float(problem.wavenumber) == np.inf:
        return np.nan * np.ones((*output_shape, 3))

    x, y, z = points.T
    k = problem.wavenumber
    h = problem.water_depth
    beta = problem.encounter_wave_direction

    wbar = x * np.cos(beta) + y * np.sin(beta)

    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        sih = np.exp(k*z)

    v = problem.g*k/problem.omega * \
        np.exp(1j * k * wbar) * \
        np.array([np.cos(problem.wave_direction) * cih, np.sin(problem.wave_direction) * cih, -1j * sih])

    return v.T.reshape((*output_shape, 3))


def airy_waves_pressure(points, problem):
    return 1j * float(problem.omega) * problem.rho * airy_waves_potential(points, problem)


def froude_krylov_force(problem):
    return problem.body.integrate_pressure(airy_waves_pressure(problem.body.mesh.faces_centers, problem))


def airy_waves_free_surface_elevation(points, problem):
    """Compute the free surface elevation at points of the undisturbed Airy waves

    Parameters
    ----------
    points: array of shape (3) or (N × 3) or (2) or (N × 2)
        coordinates of the points in which to evaluate the potential.
        If only two coordinates are passed, the last one is filled with zeros.
    problem: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    complex-valued array of shape (1,) or (N,)
        the free surface elevations
    """
    points, output_shape = _normalize_free_surface_points(points)
    return 1j * problem.omega / problem.g * airy_waves_potential(points, problem).reshape(output_shape)
