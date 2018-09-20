#!/usr/bin/env python
# coding: utf-8

import numpy as np

from capytaine.bodies import FloatingBody
from capytaine.mesh.mesh import Mesh
from capytaine.tools.geometry import Axis, Plane
from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.cylinder import HorizontalCylinder


def test_dof():
    nodes = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]])
    faces = np.array([[0, 1, 2, 3]])
    body = FloatingBody(Mesh(nodes, faces), name="one_face")
    assert body.dofs == {}

    body.add_translation_dof(direction=(1.0, 0.0, 0.0), name="1")
    assert np.allclose(body.dofs["1"], np.array([1.0, 0.0, 0.0]))

    body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    assert np.allclose(body.dofs["2"], np.array([0.0, 1.0, 0.0]))

    body.add_rotation_dof(Axis(vector=(0.0, 0.0, 1.0)), name="3")
    body.add_rotation_dof(Axis(point=(0.5, 0, 0), vector=(0.0, 0.0, 1.0)), name="4")


def test_dof_name_inference():
    body = HorizontalCylinder()
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_1'])

    body.add_rotation_dof(name="Pitch")
    body.add_rotation_dof(name="yaw")

    body.dofs.clear()
    body.add_all_rigid_body_dofs()


def test_bodies():
    body = Sphere(name="sphere", clever=False)
    assert str(body) == "sphere"
    repr(body)
    assert np.allclose(body.center, (0, 0, 0))
    body.add_translation_dof(name="Heave")

    body.extract_faces(np.where(body.mesh.faces_centers[:, 2] < 0)[0])
    body.keep_immersed_part(inplace=False)

    body.mirrored(Plane(point=(0, 1, 0), normal=(1, 0, 0)))
    body.rotated(Axis(point=(0, 1, 0), vector=(0, 0, 1)), np.pi)

    copy_of_body = body.copy(name="copy_of_sphere")
    copy_of_body.translate_x(10.0)
    copy_of_body.add_translation_dof(name="Heave")

    both = body.join_bodies(copy_of_body)
    assert set(both.dofs) == {'sphere__Heave', 'copy_of_sphere__Heave'}



