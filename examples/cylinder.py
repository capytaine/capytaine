#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation: added mass and damping of an horizontal cylinder.
"""

import numpy as np
import matplotlib.pyplot as plt

from capytaine.reference_bodies import HorizontalCylinder
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh

rho = 1000

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=80, nr=2, ntheta=20)
cylinder.translate_z(-2.0)
cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 40)

problems = [RadiationProblem(body=cylinder, rho=rho, omega=omega) for omega in omega_range]

# results = solver.solve_all(problems, processes=4)
results = [solver.solve(pb) for pb in problems]

# results = np.array(results)
# np.savetxt("results.csv", results)

# plt.figure()
# plt.plot(omega_range, results[:, 0]/(rho*cylinder.volume), label="Added mass")
# plt.plot(omega_range, results[:, 1]/(rho*cylinder.volume), label="Added damping")
# plt.legend()
# plt.show()
