#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.mesh.mesh import Mesh

from capytaine.geometric_bodies import HorizontalCylinder
from capytaine.tools.geometry import Plane

test_mesh = Mesh(vertices=np.random.rand(4, 3),faces=[range(4)], name="test_mesh")
cylinder = HorizontalCylinder().mesh.merge()

def test_init_mesh():
    dummy_mesh = Mesh()
    other_dummy_mesh = Mesh()

    assert dummy_mesh.name[:5] == "mesh_"
    assert other_dummy_mesh.name[:5] == "mesh_"
    assert int(dummy_mesh.name[5:]) + 1 == int(other_dummy_mesh.name[5:])

    with pytest.raises(AssertionError):
        Mesh(vertices=np.random.rand(6, 3), faces=[(0, 1, 2), (3, 4, 5)])

    with pytest.raises(AssertionError):
        Mesh(vertices=np.random.rand(4, 3),faces=[(0, 1, 2, np.pi)])

    with pytest.raises(AssertionError):
        Mesh(vertices=np.random.rand(3, 3),faces=[(0, 1, 2, 3)])

    assert str(test_mesh) == "Mesh(nb_vertices=4, nb_faces=1, name=test_mesh)"

def test_set_of_faces():
    assert cylinder == cylinder
    assert Mesh.from_set_of_faces(cylinder.as_set_of_faces()) == cylinder

def test_copy():
    assert test_mesh.copy() is not test_mesh
    assert test_mesh.copy() == test_mesh
    assert test_mesh.merge() is test_mesh

def test_faces():
    faces_tmp = cylinder.faces
    cylinder.faces = faces_tmp

def test_vertices():
    vertices_tmp = cylinder.vertices
    cylinder.vertices = vertices_tmp

def test_bbox():
    assert (-5, 5, -5, 5, -5, 5) == cylinder.squared_axis_aligned_bbox

def test_rotate():
    cylinder.rotate_x(np.pi)
    cylinder.rotate_y(np.pi)
    cylinder.rotate_z(np.pi)

def test_translate():
    cylinder.translate_x(0)
    cylinder.translate_y(0)
    cylinder.translate_z(0)
    cylinder.translate([0, 0, 0])

def test_merge_duplicate():
    cylinder.merge_duplicates(atol=1e-5)

def test_triangulate_quadrangles():
    cylinder.triangulate_quadrangles()

def test_mirror():
    new_cylinder = cylinder.mirror(Plane(), inplace=False)
    cylinder.mirror(Plane())
    assert new_cylinder == cylinder

def test_heal_mesh():
    cylinder.heal_mesh()

