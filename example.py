#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation.
"""

from numpy import linspace, array, savetxt

from bodies import HorizontalCylinder
from capytaine import RadiationProblem
from Nemoh import Nemoh

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=20, nr=2, ntheta=10)
cylinder.translate_z(-2.0)
cylinder.dof["Heave"] = cylinder.faces_normals @ (0, 0, 1)

solver = Nemoh()

omega_range = linspace(0.1, 5.0, 1)

problems = [RadiationProblem(bodies=[cylinder], omega=omega) for omega in omega_range]

results = []
for problem in problems:
    results.append(solver.solve(problem))

results = array(results)
savetxt("results.csv", results)

# plt.figure()
# plt.plot(omega_range, results[:, 0], label="Added mass")
# plt.plot(omega_range, results[:, 1], label="Added damping")
# plt.show()
