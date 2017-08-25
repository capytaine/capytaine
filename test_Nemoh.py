#!/usr/bin/env python
# coding: utf-8

import numpy as np
from bodies import *
from problems import RadiationProblem
from Nemoh import *

def test_floatting_sphere():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], omega=1.0, depth=np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 1819.6, rtol=1e-4)
    assert np.isclose(damping, 379.39, rtol=1e-4)

def test_floatting_sphere_finite_depth():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], omega=1.0, depth=10.0)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 1740.6, rtol=1e-4)
    assert np.isclose(damping, 380.46, rtol=1e-4)

def test_alien_sphere():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], rho=450.0, g=1.625, omega=1.0, depth=np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 515, rtol=1e-3)
    assert np.isclose(damping, 309, rtol=1e-3)

T_range, mu, nu = np.loadtxt("pytest/data/mathematica_mu_nu.tsv").T
import matplotlib.pyplot as plt
plt.figure()
plt.plot(T_range, mu, linestyle="--", label="Reference added mass")
# plt.plot(T_range, nu, linestyle="--", label="Reference added damping")

resolutions = [2, 4, 6]
for i, resolution in enumerate(resolutions):
    depth = 10.9
    flap = OpenRectangularParallelepiped(
        height=depth,
        width=3.0,
        thickness=0.01,
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
    results = np.asarray(solver.solve_all(problems, processes=4))

    plt.plot(T_range, results[:, 0], color=f'{1-(i+1)/len(resolutions)}', label=f"Added mass ({30*resolution**2} cells)")
    # plt.plot(T_range, results[:, 1], color=f'{1-(i+1)/len(resolutions)}', label="Added damping")

plt.legend()
plt.tight_layout()
plt.show()
