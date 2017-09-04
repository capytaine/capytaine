#!/usr/bin/env python
# coding: utf-8
"""
Definition of the problems to solve with the BEM.
"""

import numpy as np
from warnings import warn

from capytaine._Wavenumber import invert_xtanhx 

class RadiationProblem:
    """A radiation problem to be solved by the BEM solver."""

    def __init__(self, bodies, free_surface=0.0, sea_bottom=-np.infty, omega=1.0, rho=1000.0, g=9.81):
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

        # Clip bodies mesh
        for body in bodies:
            if any(body.vertices[:, 2] > free_surface) or any(body.vertices[:, 2] < sea_bottom):
                warn(f"""The mesh of the body {body.name} is not inside the domain.\nUse body.get_immersed_part() to clip the mesh.""")
        self.bodies = bodies

    @property
    def depth(self):
        return self.free_surface - self.sea_bottom
