#!/usr/bin/env python
# coding: utf-8
"""
Tests for the computation of the Green function and the resolution of the BEM problem.
"""

import pytest

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.mesh.meshes_collection import CollectionOfMeshes
from capytaine.mesh.symmetries import AxialSymmetry, ReflectionSymmetry, TranslationalSymmetry

from capytaine.bodies import FloatingBody
from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.cylinder import HorizontalCylinder

from capytaine.problems import RadiationProblem
from capytaine.Nemoh import Nemoh

from capytaine.tools.geometry import xOz_Plane, yOz_Plane

@pytest.mark.parametrize("reso", range(1, 3))
@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_floating_sphere(reso, depth):
    full_sphere = Sphere(radius=1.0, ntheta=reso, nphi=4*reso, clever=False, clip_free_surface=True)
    full_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=full_sphere, omega=1.0, sea_bottom=-depth)
    result1 = Nemoh().solve(problem)

    half_sphere_mesh = full_sphere.mesh.extract_faces(np.where(full_sphere.mesh.faces_centers[:, 1] > 0)[0])
    two_halves_sphere = FloatingBody(ReflectionSymmetry(half_sphere_mesh, xOz_Plane))
    two_halves_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=two_halves_sphere, omega=1.0, sea_bottom=-depth)
    result2 = Nemoh().solve(problem)

    quarter_sphere = half_sphere_mesh.extract_faces(np.where(half_sphere_mesh.faces_centers[:, 0] > 0)[0])
    four_quarter_sphere = FloatingBody(ReflectionSymmetry(ReflectionSymmetry(quarter_sphere, yOz_Plane), xOz_Plane))
    four_quarter_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=four_quarter_sphere, omega=1.0, sea_bottom=-depth)
    result3 = Nemoh().solve(problem)

    clever_sphere = Sphere(radius=1.0, ntheta=reso, nphi=4*reso, clever=True, clip_free_surface=True)
    clever_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=clever_sphere, omega=1.0, sea_bottom=-depth)
    result4 = Nemoh().solve(problem)

    # (quarter_sphere + half_sphere + full_sphere + clever_sphere).show()

    volume = 4/3*np.pi
    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result3.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result4.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result3.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result4.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)


def test_odd_axial_symmetry():
    """Buoy with odd number of slices."""
    def shape(z):
            return 0.1*(-(z+1)**2 + 16)
    buoy = FloatingBody(AxialSymmetry.from_profile(shape, z_range=np.linspace(-5.0, 0.0, 9), nphi=5))
    buoy.add_translation_dof(direction=(0, 0, 1), name="Heave")

    problem = RadiationProblem(body=buoy, omega=2.0)
    result1 = Nemoh().solve(problem)

    full_buoy = FloatingBody(buoy.mesh.merge())
    full_buoy.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=full_buoy, omega=2.0)
    result2= Nemoh().solve(problem)

    volume = buoy.mesh.volume
    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)


def test_join_translational_cylinders():
    mesh1 = HorizontalCylinder(length=10.0, radius=1.0, center=(0, 5, -5), clever=True, nr=0, ntheta=10, nx=10).mesh
    mesh2 = HorizontalCylinder(length=10.0, radius=2.0, center=(0, -5, -5), clever=True, nr=0, ntheta=10, nx=10).mesh
    joined = mesh1.join(mesh2)
    assert isinstance(joined, TranslationalSymmetry)


@pytest.mark.parametrize("depth", [10.0, np.infty])
def test_horizontal_cylinder(depth):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, clever=False, nr=2, ntheta=10, nx=10)
    assert isinstance(cylinder.mesh, Mesh)
    cylinder.translate_z(-3.0)
    cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=cylinder, omega=1.0, sea_bottom=-depth)
    result1 = Nemoh().solve(problem)

    sym_cylinder = HorizontalCylinder(length=10.0, radius=1.0, clever=True, nr=2, ntheta=10, nx=10)
    assert isinstance(sym_cylinder.mesh, CollectionOfMeshes)
    assert isinstance(sym_cylinder.mesh[0], TranslationalSymmetry)
    sym_cylinder.translate_z(-3.0)
    sym_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=sym_cylinder, omega=1.0, sea_bottom=-depth)
    result2 = Nemoh().solve(problem)

    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*cylinder.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*cylinder.volume*problem.rho)

