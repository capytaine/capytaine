#!/usr/bin/env python
# coding: utf-8
"""
Example computation: added mass and damping of an axi-symmetric buoy.
"""

import logging

import numpy as np
import matplotlib.pyplot as plt

from capytaine.reference_bodies import generate_axi_symmetric_body
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh

logging.basicConfig(level=logging.INFO, format="%(levelname)s:\t%(message)s")

rho = 1000

buoy = generate_axi_symmetric_body(
    profile=[[0, 0, -5], [1, 0, -4], [1.5, 0, -3], [2.0, 0, -2], [1.3, 0, -1], [0, 0, -0.5]]
)
# buoy.show()
buoy.dofs["Heave"] = buoy.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 40)

problems = [RadiationProblem(body=buoy, rho=rho, omega=omega) for omega in omega_range]

results = [solver.solve(pb) for pb in problems]
results = np.array(results)

plt.figure()
plt.plot(omega_range, results[:, 0, 0, 0]/(rho*buoy.volume), label="Added mass")
plt.plot(omega_range, results[:, 1, 0, 0]/(rho*buoy.volume), label="Added damping")
plt.legend()
plt.show()