#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

import numpy as np

from capytaine.reference_bodies import *
from capytaine.symmetries import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh

def test_immersed_sphere():
    sphere = Sphere(radius=1.0, ntheta=21, nphi=21)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sphere, free_surface=np.infty, sea_bottom=-np.infty)
    mass, damping = Nemoh().solve(problem)
    assert np.isclose(mass,    2187, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(damping, 0.0,  atol=1e-3*sphere.volume*problem.rho)

def test_floatting_sphere_finite_freq():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sphere, omega=1.0, sea_bottom=-np.infty)
    mass, damping = Nemoh().solve(problem)
    assert np.isclose(mass,    1819.6, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(damping, 379.39, atol=1e-3*sphere.volume*problem.rho)

def test_floatting_sphere_finite_freq_symmetry():
    half_sphere = HalfSphere(radius=1.0, ntheta=4, nphi=11, clip_free_surface=True)
    sphere = PlanarSymmetry(half_sphere, xOz_Plane)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sphere, omega=1.0, sea_bottom=-np.infty)
    mass, damping = Nemoh().solve(problem)
    assert np.isclose(mass,    1819.6, atol=1e-3*2*half_sphere.volume*problem.rho)
    assert np.isclose(damping, 379.39, atol=1e-3*2*half_sphere.volume*problem.rho)

def test_alien_sphere():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sphere, rho=450.0, g=1.625, omega=1.0, sea_bottom=-np.infty)
    mass, damping = Nemoh().solve(problem)
    assert np.isclose(mass,    515, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(damping, 309, atol=1e-3*sphere.volume*problem.rho)

def test_floatting_sphere_finite_depth():
    sphere = Sphere(radius=1.0, ntheta=7, nphi=7, clip_free_surface=True)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sphere, omega=1.0, sea_bottom=-10.0)
    mass, damping = Nemoh().solve(problem)
    assert np.isclose(mass,    1740.6, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(damping, 380.46, rtol=1e-3*sphere.volume*problem.rho)

def test_multibody():
    sphere = Sphere(radius=1.0, ntheta=11, nphi=11)
    sphere.translate_z(-2.0)
    sphere.dofs["Surge"] = sphere.faces_normals @ (1, 0, 0)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)

    cylinder = HorizontalCylinder(length=5.0, radius=1.0, nx=11, nr=2, ntheta=11)
    cylinder.translate([-1.0, 3.0, -3.0])
    cylinder.dofs["Surge"] = cylinder.faces_normals @ (1, 0, 0)
    cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

    both = cylinder + sphere
    # both.show_matplotlib()

    problem = RadiationProblem(body=both, omega=1.0, free_surface=0.0, sea_bottom=-np.infty) 
    mass, damping = Nemoh().solve(problem)

    Nemoh_2 = np.array([
        [3961.86548, 50.0367661, -3.32347107, 6.36901855E-02, 172.704819, 19.2018471, -5.67303181, -2.98873377],
        [-3.08301544, 5.72392941E-02, 14522.1689, 271.796814, 128.413834, 6.03351116, 427.167358, 64.1587067],
       [161.125534, 17.8332844, 126.392113, 5.88006783, 2242.47412, 7.17850924, 1.29002571, 0.393169671], 
       [-5.02560759, -2.75930357, 419.927460, 63.3179016, 1.23501396, 0.416424811, 2341.57593, 15.8266096], 
       ])

    assert np.allclose(mass,    Nemoh_2[:, ::2],  atol=1e-3*both.volume*problem.rho)
    assert np.allclose(damping, Nemoh_2[:, 1::2], atol=1e-3*both.volume*problem.rho)
