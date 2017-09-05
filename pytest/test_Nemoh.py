#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

import numpy as np

from capytaine.reference_bodies import *
from capytaine.symmetries import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import *

def test_panels():
    panel = OneSidedRectangle(height=1.0, width=1.0, nh=2, nw=2)
    panel.translate_z(-5.0)
    panel.rotate_x(np.pi/4)
    panel.rotate_y(np.pi/4)
    problem = RadiationProblem(bodies=[panel], omega=1.0, free_surface=0.0, sea_bottom=-np.infty)
    solver = Nemoh()
    solver.solve(problem)
    S, V = panel.build_matrices(panel)
    print(S)

test_panels()

def test_immersed_sphere():
    sphere = Sphere(radius=1.0, ntheta=21, nphi=21)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], free_surface=np.infty, sea_bottom=-np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 2187, rtol=1e-3)
    assert np.isclose(damping, 0.0, atol=1e-4)

def test_floatting_sphere_finite_freq():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], omega=1.0, sea_bottom=-np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 1819.6, rtol=1e-4)
    assert np.isclose(damping, 379.39, rtol=1e-4)

def test_floatting_sphere_finite_freq_symmetry():
    half_sphere = HalfSphere(radius=1.0, ntheta=4, nphi=11, clip_free_surface=True)
    sphere = PlanarSymmetry(half_sphere, xOz_Plane)
    # sphere.dof["Bulge"] = np.ones(sphere.nb_faces)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], omega=1.0, sea_bottom=-np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 1819.6, rtol=1e-4)
    assert np.isclose(damping, 379.39, rtol=1e-4)

#     S, V = sphere.build_matrices(sphere)
#     import matplotlib.pyplot as plt
#     plt.figure()
#     plt.imshow(np.abs(S))
#     plt.colorbar()

#     full_sphere = Sphere(radius=1.0, ntheta=7, nphi=11, clip_free_surface=True)
#     # full_sphere.dof["Bulge"] = np.ones(full_sphere.nb_faces)
#     full_sphere.dof["Heave"] = full_sphere.faces_normals @ (0, 0, 1)
#     problem = RadiationProblem(bodies=[full_sphere], omega=1.0, sea_bottom=-np.infty)
#     massf, dampingf = solver.solve(problem)
#     Sf, Vf = full_sphere.build_matrices(full_sphere)

#     plt.figure()
#     plt.imshow(np.abs(Sf))
#     plt.colorbar()
#     plt.show()

# test_floatting_sphere_finite_freq_symmetry()

def test_alien_sphere():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], rho=450.0, g=1.625, omega=1.0, sea_bottom=-np.infty)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 515, rtol=1e-3)
    assert np.isclose(damping, 309, rtol=1e-3)

def test_floatting_sphere_finite_depth():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dof["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(bodies=[sphere], omega=1.0, sea_bottom=-10.0)
    solver = Nemoh()
    mass, damping = solver.solve(problem)
    assert np.isclose(mass, 1740.6, rtol=1e-4)
    assert np.isclose(damping, 380.46, rtol=1e-4)

