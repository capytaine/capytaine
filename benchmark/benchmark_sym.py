#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation: added mass and damping of an horizontal cylinder.
"""

import cProfile, pstats, io

import numpy as np
import matplotlib.pyplot as plt

from capytaine.reference_bodies import HorizontalCylinder
from capytaine.problems import RadiationProblem
from capytaine.symmetries import *
from capytaine.Nemoh import Nemoh

omega_range = np.linspace(0.1, 5.0, 40)
rho=1000

############
#  Direct  #
############

pr = cProfile.Profile()

cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=51, nr=2, ntheta=21)
cylinder.translate_x(-5.0)
cylinder.translate_z(-2.0)
cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

problems = [RadiationProblem(body=cylinder, rho=rho, omega=omega) for omega in omega_range]
solver = Nemoh()

pr.enable()
results = [solver.solve(pb) for pb in problems]
pr.disable()

results = np.array(results)
np.savetxt("results.csv", results)

s = io.StringIO()
sortby = 'time'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
with open('c.log', 'w') as logfile:
	logfile.write(s.getvalue())

#########
#  Sym  #
#########

pr = cProfile.Profile()

quarter_cylinder = cylinder.extract_faces(
        np.where(
            np.logical_and(
                cylinder.faces_centers[:, 0] > 0,
                cylinder.faces_centers[:, 1] > 0
                )
        )[0]
        )
cylinder = PlanarSymmetry(
    PlanarSymmetry(
        quarter_cylinder,
        xOz_Plane
    ),
    yOz_Plane
)
cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

problems = [RadiationProblem(body=cylinder, rho=rho, omega=omega) for omega in omega_range]
solver = Nemoh()

pr.enable()
results = [solver.solve(pb) for pb in problems]
pr.disable()

results = np.array(results)
np.savetxt("sym_results.csv", results)

s = io.StringIO()
sortby = 'time'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
with open('symc.log', 'w') as logfile:
	logfile.write(s.getvalue())
