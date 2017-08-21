#!/usr/bin/env python
# coding: utf-8
"""
Definition of the problems to solve with the BEM.
"""

import numpy as np

from NemohCore._Wavenumber import invert_xtanhx


class RadiationProblem:
    def __init__(self, bodies, omega=1.0, depth=np.infty, rho=1000.0, g=9.81):
        self.bodies = bodies
        self.rho = rho
        self.g = g
        self.omega = omega
        self.depth = depth

        if depth == np.infty or omega**2*depth/g > 20:
            self.wavenumber = omega**2/g
        else:
            self.wavenumber = invert_xtanhx(omega**2*depth/g)/depth
