#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.bodies import FloatingBody
from capytaine.bodies_collection import CollectionOfFloatingBodies
from capytaine.symmetries import ReflectionSymmetry, xOz_Plane
from capytaine.reference_bodies import *

def test_collection():
    body_1 = generate_sphere()
    body_1.name = 'body_1'
    body_1.dofs['Dummy'] = body_1.faces_normals @ (0, 0, 0)
    body_1.dofs['Heave'] = body_1.faces_normals @ (0, 0, 1)
    assert isinstance(body_1, FloatingBody)
    assert isinstance(body_1.as_FloatingBody(), FloatingBody)

    body_2 = generate_sphere(z0=-5.0)
    body_2.name = 'body_2'
    body_2.dofs['Surge'] = body_2.faces_normals @ (1, 0, 0)
    coll = body_1 + body_2
    assert isinstance(coll, CollectionOfFloatingBodies)
    assert isinstance(coll.as_FloatingBody(), FloatingBody)
    assert coll.name == 'union_of_body_1_and_body_2'

    # Test dofs
    assert len(coll.dofs) == 3
    for name, dof in coll.dofs.items():
        assert name in ['body_1_Dummy', 'body_1_Heave', 'body_2_Surge']
        assert dof.shape[0] == coll.faces_normals.shape[0]
        assert np.all(dof[:dof.shape[0]//2] == np.zeros(coll.faces_normals.shape[0]//2)) \
                or np.all(dof[dof.shape[0]//2:] == np.zeros(coll.faces_normals.shape[0]//2))

    body_3 = generate_sphere(z0=-10.0)
    body_3.name = 'body_3'
    coll = body_1 + body_2 + body_3
    assert isinstance(coll, CollectionOfFloatingBodies)
    assert coll.name == 'union_of_body_1_body_2_and_body_3'
    assert body_1 in coll.subbodies
    assert body_2 in coll.subbodies
    assert body_3 in coll.subbodies

    assert coll.nb_vertices == coll.as_FloatingBody().nb_vertices == generate_sphere().nb_vertices*3
    assert coll.nb_faces == coll.as_FloatingBody().nb_faces == generate_sphere().nb_faces*3

    # Test the merging of identical vertices
    assert (generate_sphere() + generate_sphere()).as_FloatingBody().nb_vertices == generate_sphere().nb_vertices


def test_symmetric_bodies():
    half_sphere = generate_half_sphere(ntheta=5)
    half_sphere.name = 'half_sphere'
    full_sphere = ReflectionSymmetry(half_sphere, xOz_Plane)
    assert isinstance(full_sphere, CollectionOfFloatingBodies)
    assert half_sphere in full_sphere.subbodies
    assert full_sphere.as_FloatingBody().nb_vertices == generate_sphere(ntheta=10).nb_vertices

    other_sphere = generate_sphere(z0=-5.0)
    coll = full_sphere + other_sphere
    assert full_sphere in coll.subbodies


def test_reference_bodies():
    sphere = generate_sphere()
    # sphere.show()
    half_sphere = generate_half_sphere()
    # half_sphere.show()
    cylinder = generate_horizontal_cylinder()
    # cylinder.show()
    cylinder = generate_clever_horizontal_cylinder()
    # cylinder.show()
    rectangle = generate_one_sided_rectangle()
    # rectangle.show()
    parallelepiped = generate_open_rectangular_parallelepiped()
    # parallelepiped.show()
    parallelepiped = generate_rectangular_parallelepiped()
    # parallelepiped.show()
