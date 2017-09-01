#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from capytaine.bodies import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import *

T_range, mu, nu = np.loadtxt("pytest/data/mathematica_mu_nu.tsv").T
# plt.figure()
# plt.plot(T_range, mu, linestyle="--", label="Reference added mass")
# plt.plot(T_range, nu, linestyle="--", label="Reference added damping")

resolutions = [6, 8]
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

    # plt.plot(
    #     T_range,
    #     results[:, 0],
    #     color=f'{1-(i+1)/len(resolutions)}',
    #     label=f"Added mass ({30*resolution**2} cells)"
    # )

# plt.legend()
# plt.tight_layout()
# plt.show()
