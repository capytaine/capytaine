#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

import pytest

import numpy as np

from capytaine.geometric_bodies.sphere import generate_sphere, generate_half_sphere, generate_clever_sphere
from capytaine.reference_bodies import *
from capytaine.symmetries import *
from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh


@pytest.mark.parametrize("reso", range(1, 3))
@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_floating_sphere(reso, depth):
    full_sphere = generate_sphere(radius=1.0, ntheta=reso, nphi=4*reso, clip_free_surface=True)
    full_sphere.dofs["Heave"] = full_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=full_sphere, omega=1.0, sea_bottom=-depth)
    result1 = Nemoh().solve(problem)

    half_sphere = generate_half_sphere(radius=1.0, ntheta=reso, nphi=2*reso, clip_free_surface=True)
    # half_sphere = full_sphere.extract_faces(np.where(full_sphere.faces_centers[:, 1] > 0)[0])
    two_halves_sphere = ReflectionSymmetry(half_sphere, xOz_Plane)
    two_halves_sphere.dofs["Heave"] = two_halves_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=two_halves_sphere, omega=1.0, sea_bottom=-depth)
    result2 = Nemoh().solve(problem)

    quarter_sphere = half_sphere.extract_faces(np.where(half_sphere.faces_centers[:, 0] > 0)[0])
    four_quarter_sphere = ReflectionSymmetry(ReflectionSymmetry(quarter_sphere, yOz_Plane), xOz_Plane)
    four_quarter_sphere.dofs["Heave"] = four_quarter_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=four_quarter_sphere, omega=1.0, sea_bottom=-depth)
    result3 = Nemoh().solve(problem)

    clever_sphere = generate_clever_sphere(radius=1.0, ntheta=reso, nphi=4*reso, clip_free_surface=True)
    clever_sphere.dofs['Heave'] = clever_sphere.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=clever_sphere, omega=1.0, sea_bottom=-depth)
    result4 = Nemoh().solve(problem)

    # (quarter_sphere + half_sphere + full_sphere + clever_sphere).show()

    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result3.added_masses["Heave"], atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result4.added_masses["Heave"], atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result3.radiation_dampings["Heave"], atol=1e-4*full_sphere.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result4.radiation_dampings["Heave"], atol=1e-4*full_sphere.volume*problem.rho)


def test_odd_axial_symmetry():
    """Buoy with odd number of slices."""
    def shape(z):
            return 0.1*(-(z+1)**2 + 16)
    buoy = AxialSymmetry.from_profile(shape, z_range=np.linspace(-5.0, 0.0, 9), nphi=5)
    buoy.dofs['Heave'] = buoy.faces_normals @ (0, 0, 1)

    problem = RadiationProblem(body=buoy, omega=2.0)
    result1 = Nemoh().solve(problem)

    problem = RadiationProblem(body=buoy.as_FloatingBody(), omega=2.0)
    result2= Nemoh().solve(problem)

    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*buoy.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*buoy.volume*problem.rho)


@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_horizontal_cylinder(depth):
    cylinder = generate_open_horizontal_cylinder(length=10.0, radius=1.0, ntheta=10, nx=10)
    cylinder.translate_z(-3.0)
    cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=cylinder, omega=1.0, sea_bottom=-depth)
    result1 = Nemoh().solve(problem)

    sym_cylinder = generate_clever_horizontal_cylinder(length=10.0, radius=1.0, ntheta=10, nx=10)
    sym_cylinder.translate_z(-3.0)
    sym_cylinder.dofs["Heave"] = sym_cylinder.faces_normals @ (0, 0, 1)
    problem = RadiationProblem(body=sym_cylinder, omega=1.0, sea_bottom=-depth)
    result2 = Nemoh().solve(problem)

    cylinder_volume = 10*1.0*2*np.pi
    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*cylinder_volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*cylinder_volume*problem.rho)

