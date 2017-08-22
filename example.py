#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation.
"""

import numpy as np
import matplotlib.pyplot as plt

from bodies import HorizontalCylinder
from capytaine import RadiationProblem
from Nemoh import Nemoh

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=20, nr=2, ntheta=10)
cylinder.translate_z(-2.0)
cylinder.dof["Heave"] = cylinder.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = np.linspace(0.1, 5.0, 10)

problems = [RadiationProblem(bodies=[cylinder], omega=omega) for omega in omega_range]

results = []
for problem in problems:
    results.append(solver.solve(problem))

results = np.array(results)

plt.figure()
plt.plot(omega_range, results[:, 0], label="Added mass")
plt.plot(omega_range, results[:, 1], label="Added damping")
plt.legend()
plt.show()
