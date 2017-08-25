#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation.
"""

import numpy as np
import matplotlib.pyplot as plt

from bodies import HorizontalCylinder
from problems import RadiationProblem
from Nemoh import Nemoh

rho = 1000

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=80, nr=2, ntheta=20)
cylinder.translate_z(-2.0)
cylinder.dof["Heave"] = cylinder.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 40)

problems = [RadiationProblem(bodies=[cylinder], rho=rho, depth=np.infty, omega=omega) for omega in omega_range]

results = solver.solve_all(problems, processes=4)

results = np.array(results)
# np.savetxt("results.csv", results)

plt.figure()
plt.plot(omega_range, results[:, 0]/(rho*cylinder.volume), label="Added mass")
plt.plot(omega_range, results[:, 1]/(rho*cylinder.volume), label="Added damping")
plt.legend()
plt.show()
