#!/usr/bin/env python
# coding: utf-8

import capytaine
from capytaine.bodies import HorizontalCylinder

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nlength=10, nradius=2)
cylinder.translate_mesh(z=-5.0)
cylinder.dof["Heave"] = np.dot(cylinder.mesh.normal_vectors, (0,0,-1))

other_body = capytaine.load_mesh("/path/to/mesh")

test_case = capytaine.RadiationProblem(
    bodies=[cylinder, otherBody],
    frequency=0.1
    rho=1000.0,
    g=9.81
)

capytaine.Nemoh().solve(test_case, save_pressure=True)

A = test_case.added_mass_matrix
