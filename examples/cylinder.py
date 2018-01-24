#!/usr/bin/env python
# coding: utf-8
"""
Example computation: added mass and damping of an horizontal cylinder.
"""

import logging

import numpy as np
import matplotlib.pyplot as plt

from capytaine.reference_bodies import generate_horizontal_cylinder
from capytaine.problems import RadiationProblem
from capytaine.results import assemble_radiation_results_matrices
from capytaine.Nemoh import Nemoh

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

rho = 1000

cylinder = generate_horizontal_cylinder(length=10.0, radius=1.0, nx=20, nr=2, ntheta=20)
cylinder.translate_z(-2.0)
cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 40)
problems = [RadiationProblem(body=cylinder, radiating_dof="Heave", rho=rho, omega=omega) for omega in omega_range]

results = [solver.solve(pb) for pb in problems]

added_masses, radiation_dampings = assemble_radiation_results_matrices(results)
np.savetxt("added_masses.csv", added_masses[:, 0, 0])

plt.figure()
plt.plot(omega_range, added_masses[:, 0, 0]/(rho*cylinder.volume), label="Added mass")
plt.plot(omega_range, radiation_dampings[:, 0, 0]/(rho*cylinder.volume), label="Added damping")
plt.legend()
plt.show()
