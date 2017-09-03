#!/usr/bin/env python
# coding: utf-8
"""
Exemple computation: added mass and damping of an oscillating flap.
"""

import numpy as np
import matplotlib.pyplot as plt

from capytaine.bodies import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import *

T_range, mu, nu = np.loadtxt("data/flap_mu_nu.tsv").T

resolutions = [2]
for i, resolution in enumerate(resolutions):
    depth = 10.9
    flap = OpenRectangularParallelepiped(
        height=depth,
        width=3.0,
        thickness=0.001,
        nh=int(3*resolution),
        nw=int(10*resolution),
        nth=2
    )
    flap.translate_z(-depth)
    flap.dof["Oscillation"] = np.asarray([
        flap.faces_normals[j, 1] *
        (flap.faces_centers[j, 2] + 9.4) * np.heaviside(flap.faces_centers[j, 2] + 9.4, 0.0)
        for j in range(flap.nb_faces)])

    problems = [RadiationProblem(bodies=[flap], omega=omega, depth=depth) for omega in 2*np.pi/T_range]
    solver = Nemoh()
    results = np.asarray(solver.solve_all(problems, processes=3))

    np.savetxt(f"Results_{30*resolution**2}_cells.tsv", np.asarray(results))

