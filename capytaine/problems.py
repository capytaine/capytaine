#!/usr/bin/env python
# coding: utf-8
"""
Definition of the problems to solve with the BEM.
"""

from warnings import warn

import numpy as np

from capytaine._Wavenumber import invert_xtanhx


class PotentialFlowProblem:

    def __init__(self, body, free_surface=0.0, sea_bottom=-np.infty, omega=1.0, rho=1000.0, g=9.81):
        self.rho = rho
        self.g = g
        self.omega = omega

        if free_surface < sea_bottom:
            raise Exception("Sea bottom is above the free surface.")

        self.free_surface = free_surface
        self.sea_bottom = sea_bottom

        if self.depth == np.infty or omega**2*self.depth/g > 20:
            self.wavenumber = omega**2/g
        else:
            self.wavenumber = invert_xtanhx(omega**2*self.depth/g)/self.depth

        if any(body.vertices[:, 2] > free_surface) or any(body.vertices[:, 2] < sea_bottom):
            warn(f"""The mesh of the body {body.name} is not inside the domain.\nUse body.get_immersed_part() to clip the mesh.""")
        self.body = body

    @property
    def depth(self):
        return self.free_surface - self.sea_bottom

    @property
    def wavelength(self):
        return 2*np.pi/self.wavenumber

    @property
    def period(self):
        return 2*np.pi/self.omega


class DiffractionProblem(PotentialFlowProblem):

    def __init__(self, *args, angle=0.0, **kwargs):
        self.angle = angle
        PotentialFlowProblem.__init__(self, *args, **kwargs)

    def __str__(self):
        return f"Diffraction problem of {self.body.name} with depth={self.free_surface-self.sea_bottom:.1e}, angle={self.angle:.3f} and omega={self.omega:.3f}"

    def __repr__(self):
        return f"DiffractionProblem(body={self.body.name}, free_surface={self.free_surface}, sea_bottom={self.sea_bottom}, angle={self.angle}, omega={self.omega}, rho={self.rho}, g={self.g})"

    def Airy_wave_potential(self, X):
        """Compute the potential for Airy waves at a given point (or array of points).

        Parameters
        ----------
        X: array (3) or (N x 3)
            The coordinates of the points in which to evaluate the potential.

        Returns
        -------
        array (1) or (N x 1)
            The potential
        """
        x, y, z = X.T
        k = self.wavenumber
        h = self.depth
        wbar = x*np.cos(self.angle) + y*np.sin(self.angle)

        if k*h < 20 and k*h >= 0:
            cih = np.cosh(k*(z+h))/np.cosh(k*h)
            sih = np.sinh(k*(z+h))/np.cosh(k*h)
        else:
            cih = np.exp(k*z)
            sih = np.exp(k*z)

        return -1j*self.g/self.omega * cih * np.exp(1j * k * wbar)

    def Airy_wave_velocity(self, X):
        """Compute the fluid velocity for Airy waves at a given point (or array of points).

        Parameters
        ----------
        X: array (3) or (N x 3)
            The coordinates of the points in which to evaluate the velocity.

        Returns
        -------
        array (3) or (N x 3)
            The velocity vectors
        """
        x, y, z = X.T
        k = self.wavenumber
        h = self.depth

        wbar = x*np.cos(self.angle) + y*np.sin(self.angle)

        if k*h < 20 and k*h >= 0:
            cih = np.cosh(k*(z+h))/np.cosh(k*h)
            sih = np.sinh(k*(z+h))/np.cosh(k*h)
        else:
            cih = np.exp(k*z)
            sih = np.exp(k*z)

        v = self.g*k/self.omega * \
                np.exp(1j * k * wbar) * \
                np.array([np.cos(self.angle)*cih, np.sin(self.angle)*cih, -1j*sih])

        return v.T


class RadiationProblem(PotentialFlowProblem):
    """A radiation problem to be solved by the BEM solver."""

    def __init__(self, *args, angle=0.0, **kwargs):
        self.sources = {}
        self.potential = {}
        PotentialFlowProblem.__init__(self, *args, **kwargs)

    def __str__(self):
        return f"Radiation problem of {self.body.name} with depth={self.free_surface-self.sea_bottom:.1e} and omega={self.omega:.3f}"

    def __repr__(self):
        return f"RadiationProblem(body={self.body.name}, free_surface={self.free_surface}, sea_bottom={self.sea_bottom}, omega={self.omega}, rho={self.rho}, g={self.g})"

    @property
    def dofs(self):
        return self.body.dofs

