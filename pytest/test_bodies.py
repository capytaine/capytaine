#!/usr/bin/env python
# coding: utf-8

import numpy as np

from capytaine.bodies import FloatingBody
from capytaine.mesh.mesh import Mesh
from capytaine.geometric_bodies.sphere import Sphere


def test_dof():
    nodes = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]])
    faces = np.array([[0, 1, 2, 3]])
    body = FloatingBody(Mesh(nodes, faces), name="one_face")
    assert body.dofs == {}

    # body.add_translation_dof(direction=(1.0, 0.0, 0.0), name="1")
    # assert body.dofs["1"] == np.array([1.0, 0.0, 0.0])
    #
    # body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    # assert body.dofs["2"] == np.array([0.0, 1.0, 0.0])
    #
    # body.add_rotation_dof(axis_direction=(0.0, 0.0, 1.0), name="3")
    # assert body.dofs["3"] == np.array([0.5])
    #
    # body.add_rotation_dof(axis_point=(0.5, 0, 0), axis_direction=(0.0, 0.0, 1.0), name="4")
    # assert body.dofs["4"] == np.array([0.0])


def test_dof_name_inference():
    body = Sphere()
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_1'])

