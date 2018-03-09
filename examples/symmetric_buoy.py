#!/usr/bin/env python
# coding: utf-8
"""
Example computation: added mass and damping of an axi-symmetric buoy.
"""

import logging

import numpy as np
import matplotlib.pyplot as plt

from capytaine import *

logging.basicConfig(level=logging.INFO, format="%(levelname)s:\t%(message)s")

rho = 1000

def shape(z):
    return 0.1*(-(z+1)**2 + 16)

# buoy = generate_axi_symmetric_body(
#     profile=[[0, 0, -5], [1, 0, -4], [1.5, 0, -3], [2.0, 0, -2], [1.3, 0, -1], [0, 0, -0.5]]
# )

buoy = generate_axi_symmetric_body(shape, z_range=np.linspace(-5.0, 0.0, 30), nphi=40)
# buoy.show()

buoy.dofs["Heave"] = buoy.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 60)

problems = [RadiationProblem(body=buoy, radiating_dof='Heave', rho=rho, omega=omega) for omega in omega_range]

results = [solver.solve(pb) for pb in problems]
dataset = assemble_dataset(results)

plt.figure()
plt.plot(omega_range,
         dataset['added_mass'].sel(radiating_dof='Heave', influenced_dof='Heave')/(rho*buoy.volume),
         label="Added mass")
plt.plot(omega_range,
         dataset['radiation_damping'].sel(radiating_dof='Heave', influenced_dof='Heave')/(rho*buoy.volume),
         label="Added damping")
plt.grid()
plt.legend()
plt.show()
