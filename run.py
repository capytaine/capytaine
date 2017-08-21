#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation.
"""

from bodies import HorizontalCylinder
from capytaine import RadiationProblem
from Nemoh import Nemoh

cylinder = HorizontalCylinder(length=1.0, radius=1.0, nx=10, nr=2, ntheta=10)
cylinder.translate_z(-2.0)
cylinder.dof["Heave"] = cylinder.faces_normals @ (0, 0, 1)

test_case = RadiationProblem(bodies=[cylinder], omega=0.1)

solver = Nemoh()
mass, damping = solver.solve(test_case)
print(mass, damping)
