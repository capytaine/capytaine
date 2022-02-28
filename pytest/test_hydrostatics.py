#!/usr/bin/env python
# coding: utf-8
"""
Tests for the functions that computes hydrostatic from the mesh vertices
and faces
"""

import json
from pathlib import Path

import capytaine as cpt
import numpy as np

def test_all_hydrostatics():
    density = 1000
    gravity = 9.80665

    sphere = cpt.Sphere(
        radius=10.0,          # Dimension
        center=(0,0,-1),   # Position
        nphi=100, ntheta=50,  # Fineness of the mesh
    )
    horizontal_cylinder = cpt.HorizontalCylinder(
        length=10.0, radius=5.0,  # Dimensions
        center=(0,10,-1),        # Position
        nr=10, nx=10, ntheta=10,   # Fineness of the mesh
    )
    vertical_cylinder = cpt.VerticalCylinder(
        length=10.0, radius=5.0,  # Dimensions
        center=(0,0,0),        # Position
        nr=10, nx=10, ntheta=10,   # Fineness of the mesh
    )
    cog = np.zeros(3)
    body = sphere + horizontal_cylinder + vertical_cylinder
    body.add_all_rigid_body_dofs()

    capy_hsdb = body.compute_hydrostatics(cog=cog, density=density, gravity=gravity)
    capy_hsdb["stiffness_matrix"] = capy_hsdb["stiffness_matrix"][2:5,2:5]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"][3:,3:]
    # =============================================================================
    # Meshmagick
    # =============================================================================
    case_dir = Path(__file__).parent / "Hydrostatics_cases"
    # import meshmagick.mesh as mmm
    # import meshmagick.hydrostatics as mmhs
    # body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)
    # mm_hsdb = mmhs.compute_hydrostatics(body_mesh, np.array(cog), density, gravity)
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
            assert np.isclose(capy_hsdb[var], mm_hsdb[var], 
                              rtol=1e-2, atol=1e-3).all()

def test_radial_elastic_dof():
    body = cpt.Sphere(
    radius=10.0,          # Dimension
    center=(0,0,-10),   # Position
    nphi=100, ntheta=100,  # Fineness of the mesh
    ) 
    body.dofs["radial"] = body.mesh.faces_normals
    capy_hs = body.get_hydrostatic_stiffness()
    analytical_hs = 0.0

    assert np.isclose(capy_hs, analytical_hs)

def test_vertical_elastic_dof():
    
    length = 10.0
    radius = 5.0
    center = (0.0,0.0,0.0)
    cylinder_keel = length/2 - center[2]
    
    body = cpt.VerticalCylinder(
        length=length, radius=radius,  # Dimensions
        center=center,        # Position
        nr=100, nx=100, ntheta=100,   # Fineness of the mesh
    )
    body.keep_immersed_part()
    faces_centers = body.mesh.faces_centers
    body.dofs["elongate"] = np.zeros_like(faces_centers)
    body.dofs["elongate"][:,2] = faces_centers[:,2]
    density = 1000
    gravity = 9.80665
    capy_hs = body.get_hydrostatic_stiffness(density=density, gravity=gravity)
    analytical_hs = np.pi * radius**2 * density * gravity * cylinder_keel**2

    assert np.isclose(capy_hs, analytical_hs, rtol=1e-3)


