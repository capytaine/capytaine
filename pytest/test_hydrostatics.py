#!/usr/bin/env python
# coding: utf-8
"""
Tests for the functions that computes hydrostatic from the mesh vertices 
and faces
"""

import pytest

import capytaine as cpt
import numpy as np
import json
import pprint 
from capytaine.bodies import FloatingBody
import meshmagick.mesh as mmm
import meshmagick.hydrostatics as mmhs


def test_sphere():
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

    cog = (0,0,0)
    old_body = sphere + horizontal_cylinder + vertical_cylinder
    body = FloatingBody(mesh=old_body.mesh, name="Pod")
    body.add_all_rigid_body_dofs()

    capy_hsdb = body.compute_hydrostatics(cog=np.array(cog), density=density, gravity=gravity)
    capy_hsdb["stiffness_matrix"] = capy_hsdb["stiffness_matrix"][2:5,2:5]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"][3:,3:]
    # =============================================================================
    # Meshmagick
    # =============================================================================

    body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)
    mm_hsdb = mmhs.compute_hydrostatics(body_mesh, np.array(cog), density, gravity)

    mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(rho_medium=density).inertia_matrix

    # =============================================================================
    # Logging
    # =============================================================================
    for var in capy_hsdb.keys():
        if var in mm_hsdb.keys():
            assert(np.isclose(capy_hsdb[var], mm_hsdb[var], rtol=1e-2, 
                              atol=1e-3).all())