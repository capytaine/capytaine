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

    def __init__(self, bodies, omega=1.0, depth=np.infty, rho=1000.0, g=9.81):
        self.rho = rho
        self.g = g
        self.omega = omega
        self.depth = depth

        if depth == np.infty or omega**2*depth/g > 20:
            self.wavenumber = omega**2/g
        else:
            self.wavenumber = invert_xtanhx(omega**2*depth/g)/depth

        # Clip bodies mesh
        for body in bodies:
            if any(body.vertices[:, 2] > 1e-5) or any(body.vertices[:, 2] < -depth):
                warn(f"""The mesh of the body {body.name} is not inside the domain.\nUse body.get_immerged_part() to clip the mesh.""")
        self.bodies = bodies

