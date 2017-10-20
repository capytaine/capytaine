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


class DiffractionProblem(PotentialFlowProblem):

    def __init__(self, *args, angle=0.0, **kwargs):
        self.angle = angle
        PotentialFlowProblem.__init__(self, *args, **kwargs)

    def __str__(self):
        return f"Diffraction problem of {self.body.name} with depth={self.free_surface-self.sea_bottom:.1e}, angle={self.angle:.3f} and omega={self.omega:.3f}"

    def __repr__(self):
        return f"DiffractionProblem(body={self.body.name}, free_surface={self.free_surface}, sea_bottom={self.sea_bottom}, angle={self.angle}, omega={self.omega}, rho={self.rho}, g={self.g})"

    def airy_wave(self, X):

        x, y, z = X

        XEFF, YEFF = 0, 0
        wbar = (x-XEFF)*np.cos(self.angle) + (y-YEFF)*np.sin(self.angle)

        if self.wavenumber*self.depth < 20 and self.wavenumber*self.depth >= 0:
            cih = np.cosh(self.wavenumber*(z+self.depth))/np.cosh(self.wavenumber*self.depth)
            sih = np.sinh(self.wavenumber*(z+self.depth))/np.sinh(self.wavenumber*self.depth)
        else:
            cih = np.exp(self.wavenumber*z)
            sih = np.exp(self.wavenumber*z)

        v = self.g*self.wavenumber/self.omega * \
                np.array([np.cos(self.angle)*cih, np.sin(self.angle)*cih, -1j*sih]) * \
                np.exp(1j * self.wavenumber * wbar)

        p = self.rho * self.g * cih * np.exp(1j * self.wavenumber * wbar)

        return p, v


class RadiationProblem(PotentialFlowProblem):
    """A radiation problem to be solved by the BEM solver."""

    def __str__(self):
        return f"Radiation problem of {self.body.name} with depth={self.free_surface-self.sea_bottom:.1e} and omega={self.omega:.3f}"

    def __repr__(self):
        return f"RadiationProblem(body={self.body.name}, free_surface={self.free_surface}, sea_bottom={self.sea_bottom}, omega={self.omega}, rho={self.rho}, g={self.g})"

    @property
    def dofs(self):
        return self.body.dofs

