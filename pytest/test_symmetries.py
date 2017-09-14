#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

from collections import Counter
import pytest

import numpy as np

from capytaine.reference_bodies import *
from capytaine.symmetries import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh

def test_panels():
    panel = OneSidedRectangle(height=1.0, width=1.0, nh=7, nw=3)
    panel.translate_z(-1.0)
    half_panel = panel.extract_faces(np.where(panel.faces_centers[:, 0] > 0)[0])
    symmetric_panel = ReflectionSymmetry(half_panel, yOz_Plane)

    # symmetric_panel.show_matplotlib()

    # Next lines only to set up LISC...
    problem = RadiationProblem(body=panel, omega=0.1, free_surface=0.0, sea_bottom=-np.infty)
    Nemoh().solve(problem)

    S1, V1 = panel.build_matrices(panel)
    S2, V2 = symmetric_panel.build_matrices(symmetric_panel)

    # import matplotlib.pyplot as plt
    # plt.matshow(np.real(S1), vmin=-0.1, vmax=0)
    # plt.colorbar()
    # plt.matshow(np.real(S2), vmin=-0.1, vmax=0)
    # plt.colorbar()
    # plt.show()

    assert np.allclose(S1, S2, atol=1e-5)
    assert np.allclose(V1, V2, atol=1e-5)

@pytest.mark.parametrize("reso", range(2, 5))
def test_floating_sphere(reso):
    full_sphere = Sphere(radius=1.0, ntheta=2*reso+1, nphi=2*reso+1, clip_free_surface=True)
    full_sphere.dofs["Heave"] = full_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=full_sphere, omega=1.0, sea_bottom=-np.infty)
    mass1, damping1 = Nemoh().solve(problem)

    half_sphere = HalfSphere(radius=1.0, ntheta=reso+1, nphi=2*reso+1, clip_free_surface=True)
    # half_sphere = full_sphere.extract_faces(np.where(full_sphere.faces_centers[:, 1] > 0)[0])
    two_halves_sphere = ReflectionSymmetry(half_sphere, xOz_Plane)
    two_halves_sphere.dofs["Heave"] = two_halves_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=two_halves_sphere, omega=1.0, sea_bottom=-np.infty)
    mass2, damping2 = Nemoh().solve(problem)

    quarter_sphere = half_sphere.extract_faces(np.where(half_sphere.faces_centers[:, 0] > 0)[0])
    four_quarter_sphere = ReflectionSymmetry(ReflectionSymmetry(quarter_sphere, yOz_Plane), xOz_Plane)
    four_quarter_sphere.dofs["Heave"] = four_quarter_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=four_quarter_sphere, omega=1.0, sea_bottom=-np.infty)
    mass3, damping3 = Nemoh().solve(problem)

    # (quarter_sphere + half_sphere + full_sphere).show_matplotlib()

    assert np.isclose(mass1,    mass2,    atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(damping1, damping2, atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(mass1,    mass3,    atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(damping1, damping3, atol=1e-4*full_sphere.volume*problem.rho)

def test_horizontal_cylinder():
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, ntheta=11, nr=0, nx=11)
    cylinder.translate_z(-3.0)
    cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=cylinder, omega=1.0, sea_bottom=-np.infty)
    mass1, damping1 = Nemoh().solve(problem)

    ring = HorizontalCylinder(length=1.0, radius=1.0, ntheta=11, nr=0, nx=2)
    ring.translate_z(-3.0)
    sym_cylinder = TranslationalSymmetry(ring, translation=(1.0, 0.0, 0.0), nb_repetitions=9)
    sym_cylinder.dofs["Heave"] = sym_cylinder.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sym_cylinder, omega=1.0, sea_bottom=-np.infty)
    mass2, damping2 = Nemoh().solve(problem)

    assert np.isclose(mass1,    mass2,    atol=1e-4*cylinder.volume*problem.rho)
    assert np.isclose(damping1, damping2, atol=1e-4*cylinder.volume*problem.rho)

# print(Counter(np.around(np.real(S1.flatten()), decimals=2)))
# print(Counter(np.around(np.real(S3.flatten()), decimals=2)))

