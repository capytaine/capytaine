#!/usr/bin/env python
# coding: utf-8
"""
Example computation: added mass and damping of an horizontal cylinder.
"""

import logging

import numpy as np
import matplotlib.pyplot as plt

from capytaine import *

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

rho = 1000

cylinder = generate_clever_horizontal_cylinder(length=10.0, radius=1.0, nx=30, ntheta=30)
cylinder.translate_z(-2.0)
cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)
cylinder.dofs["Sway"] = cylinder.faces_normals @ (0, 1, 0)

solver = Nemoh()

omega_range = np.linspace(0.1, 6.0, 40)
problems = [RadiationProblem(body=cylinder, radiating_dof=dof, rho=rho, omega=omega)
            for dof in cylinder.dofs for omega in omega_range]

results = [solver.solve(pb) for pb in sorted(problems)]
data = assemble_dataset(results)

# np.savetxt("added_masses.csv", data['added_mass'].sel(radiating_dof='Heave', influenced_dof='Heave'))

plt.figure()
plt.plot(omega_range,
         data['added_mass'].sel(radiating_dof='Heave', influenced_dof='Heave')/(rho*cylinder.volume),
         label="Added mass")
plt.plot(omega_range,
         data['radiation_damping'].sel(radiating_dof='Heave', influenced_dof='Heave')/(rho*cylinder.volume),
         label="Added damping")
plt.legend()
plt.show()
