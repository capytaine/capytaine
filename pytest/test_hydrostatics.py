#!/usr/bin/env python
# coding: utf-8
"""
Tests for the functions that computes hydrostatic from the mesh vertices
and faces
"""

import logging
import pytest
import json
from pathlib import Path

import capytaine as cpt
import numpy as np

# TODO: Can I use pytest fixtures to avoid regenerating so many spheres?
# Does pytest fixture do a copy of the object, since I modifying the sphere in-place?

def test_disp_mass_of_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=50, ntheta=50)
    analytical_volume = 4/3*np.pi*1.0**3
    assert np.isclose(sphere.disp_mass(rho=1000), 1000*analytical_volume, rtol=2e-2)

def test_waterplane_area_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    assert np.isclose(sphere.waterplane_area, 0.0)

def test_waterplane_center_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    assert sphere.waterplane_center is None

def test_waterplane_center_of_sphere_at_surface():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    assert np.allclose(sphere.waterplane_center, [0.0, 0.0])

def test_stiffness_when_no_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    with pytest.raises(AttributeError, match=".* no dof .*"):
        sphere.compute_hydrostatic_stiffness()

def test_stiffness_with_divergence():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.dofs["elongate_in_z"] = np.array([(0, 0, z) for (x, y, z) in sphere.mesh.faces_centers])
    hs_1 = sphere.compute_hydrostatic_stiffness()
    hs_2 = sphere.compute_hydrostatic_stiffness(divergence={"elongate_in_z": np.ones(sphere.mesh.nb_faces)})
    assert hs_1.values[0, 0] != hs_2.values[0, 0]
    analytical_hs = - 1000.0 * 9.81 * (4 * sphere.volume * sphere.center_of_buoyancy[2])
    assert np.isclose(hs_2.values[0, 0], analytical_hs)

def test_stiffness_with_malformed_divergence(caplog):
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.dofs["elongate_in_z"] = np.array([(0, 0, z) for (x, y, z) in sphere.mesh.faces_centers])
    hs_1 = sphere.compute_hydrostatic_stiffness()
    with caplog.at_level(logging.WARNING):
        hs_2 = sphere.compute_hydrostatic_stiffness(divergence={"foobar": np.ones(sphere.mesh.nb_faces)})
    assert hs_1.values[0, 0] == hs_2.values[0, 0]
    assert "without the divergence" in caplog.text

def test_mass_of_sphere_for_non_default_density():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=50, ntheta=50)
    sphere.add_translation_dof(name="Heave")
    sphere.center_of_mass = np.array([0, 0, -2])
    m = sphere.compute_rigid_body_inertia(rho=500)
    analytical_volume = 4/3*np.pi*1.0**3
    assert np.isclose(m.values[0, 0], 500*analytical_volume, rtol=2e-2)

def test_inertia_rigid_body_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.add_all_rigid_body_dofs()
    assert np.all(sphere.compute_rigid_body_inertia(output_type="rigid_dofs")
            == sphere.compute_rigid_body_inertia(output_type="all_dofs"))
    assert np.all(sphere.compute_rigid_body_inertia(output_type="rigid_dofs")
            == sphere.compute_rigid_body_inertia(output_type="body_dofs"))

def test_inertia_wrong_output_type():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.add_all_rigid_body_dofs()
    with pytest.raises(ValueError):
        sphere.compute_rigid_body_inertia(output_type="foo")

def test_inertia_when_no_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    m =sphere.compute_rigid_body_inertia()
    assert m.shape == (0, 0)

def test_hydrostatics_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    sphere.add_all_rigid_body_dofs()
    sphere.center_of_mass = np.array([0, 0, -2])
    sphere.compute_hydrostatics()


def test_all_hydrostatics():
    density = 1000
    gravity = 9.80665

    sphere = cpt.Sphere(
        radius=10.0,
        center=(0,0,-1),
        nphi=100, ntheta=50,
    )
    horizontal_cylinder = cpt.HorizontalCylinder(
        length=10.0, radius=5.0,
        center=(0,10,-1),
        nr=100, nx=100, ntheta=10,
    )
    vertical_cylinder = cpt.VerticalCylinder(
        length=10.0, radius=5.0,
        center=(10,0,0),
        nr=200, nx=200, ntheta=10,
    )
    body = sphere + horizontal_cylinder + vertical_cylinder
    body.add_all_rigid_body_dofs()
    body.center_of_mass = body.center_of_buoyancy

    capy_hsdb = body.compute_hydrostatics(rho=density, g=gravity)

    stiff_compare_dofs = ["Heave", "Roll", "Pitch"]
    capy_hsdb["stiffness_matrix"] = capy_hsdb["hydrostatic_stiffness"].sel(
        influenced_dof=stiff_compare_dofs, radiating_dof=stiff_compare_dofs
        ).values

    mass_compare_dofs = ["Roll", "Pitch", "Yaw"]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"].sel(
        influenced_dof=mass_compare_dofs, radiating_dof=mass_compare_dofs
        ).values


    # =============================================================================
    # Meshmagick
    # =============================================================================
    case_dir = Path(__file__).parent / "Hydrostatics_cases"
    # case_dir.mkdir(parents=True, exist_ok=True)
    # import meshmagick.mesh as mmm
    # import meshmagick.hydrostatics as mmhs
    # body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)
    # mm_hsdb = mmhs.compute_hydrostatics(body_mesh, body.center_of_mass, density, gravity)
    # mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(
    # rho_medium=density).inertia_matrix
    # mm_hsdb["mesh"] = ""
    # with open(f'{case_dir}/sphere__hor_cyl__ver_cyl.pkl.json', 'w') as convert_file:
    #     mm_hsdb_json = {key:(value.tolist() if type(value)==np.ndarray else value)
    #                         for key, value in mm_hsdb.items() }
    #     convert_file.write(json.dumps(mm_hsdb_json))

    with open(f'{case_dir}/sphere__hor_cyl__ver_cyl.pkl.json', 'r',
              encoding="UTF-8") as f:
        mm_hsdb = json.load(f)

    # =============================================================================
    # Testing
    # =============================================================================

    for var in capy_hsdb.keys():
        if var in mm_hsdb.keys():
            # if not np.isclose(capy_hsdb[var], mm_hsdb[var],
            #                   rtol=1e-2, atol=1e-3).all():
            #     print(f"{var}:")
            #     print(f"    Capytaine  - {capy_hsdb[var]}")
            #     print(f"    Meshmagick - {mm_hsdb[var]}")
            assert np.isclose(capy_hsdb[var], mm_hsdb[var],
                              rtol=1e-2, atol=1e-3).all()


def test_vertical_elastic_dof():

    bodies = [
        cpt.VerticalCylinder(
            length=10.0, radius=5.0,
            center=(0.0,0.0,0.0),
            nr=100, nx=100, ntheta=100,
        ),

        cpt.Sphere(
            radius=100,
            center=(0,0,0),
            nphi=50, ntheta=50,
        ),

        cpt.HorizontalCylinder(
            length=5.0, radius=1.0,
            center=(0,10,0),
            nr=20, nx=20, ntheta=10,
        )
    ]

    for body in bodies:
        body.center_of_mass = body.center_of_buoyancy
        body.keep_immersed_part()
        faces_centers = body.mesh.faces_centers
        body.dofs["elongate"] = np.zeros_like(faces_centers)
        body.dofs["elongate"][:,2] = faces_centers[:,2]

        divergence = np.ones(body.mesh.faces_centers.shape[0])

        density = 1000
        gravity = 9.80665
        capy_hs = body.each_hydrostatic_stiffness("elongate", "elongate",
                    influenced_dof_div=divergence, rho=density, g=gravity).values[0][0]

        analytical_hs = - (density * gravity * 4 * body.volume
                        * body.center_of_buoyancy[2])

        assert np.isclose(capy_hs, analytical_hs)


