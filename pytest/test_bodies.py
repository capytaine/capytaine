#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.bodies import FloatingBody
from capytaine.bodies_collection import CollectionOfFloatingBodies
from capytaine.symmetries import ReflectionSymmetry, xOz_Plane
from capytaine.reference_bodies import *

def test_collection():
    body_1 = Sphere()
    body_1.name = 'body_1'
    body_1.dofs['Dummy'] = body_1.faces_normals @ (0, 0, 0)
    body_1.dofs['Heave'] = body_1.faces_normals @ (0, 0, 1)
    assert isinstance(body_1, FloatingBody)
    assert isinstance(body_1.as_FloatingBody(), FloatingBody)

    body_2 = Sphere(z0=-5.0)
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

    body_3 = Sphere(z0=-10.0)
    body_3.name = 'body_3'
    coll = body_1 + body_2 + body_3
    assert isinstance(coll, CollectionOfFloatingBodies)
    assert coll.name == 'union_of_body_1_body_2_and_body_3'
    assert body_1 in coll.subbodies
    assert body_2 in coll.subbodies
    assert body_3 in coll.subbodies

    assert coll.nb_vertices == coll.as_FloatingBody().nb_vertices == Sphere().nb_vertices*3
    assert coll.nb_faces == coll.as_FloatingBody().nb_faces == Sphere().nb_faces*3

    # Test the merging of identical vertices
    assert (Sphere() + Sphere()).as_FloatingBody().nb_vertices == Sphere().nb_vertices


def test_symmetric_bodies():
    half_sphere = HalfSphere(ntheta=6)
    half_sphere.name = 'half_sphere'
    full_sphere = ReflectionSymmetry(half_sphere, xOz_Plane)
    assert isinstance(full_sphere, CollectionOfFloatingBodies)
    assert half_sphere in full_sphere.subbodies
    assert full_sphere.as_FloatingBody().nb_vertices == Sphere(ntheta=11).nb_vertices

    other_sphere = Sphere(z0=-5.0)
    coll = full_sphere + other_sphere
    assert full_sphere in coll.subbodies


@pytest.mark.parametrize("size", np.linspace(1, 10, 2))
@pytest.mark.parametrize("ncells", [6, 11])
def test_parallelepiped(size, ncells):
    rp = RectangularParallelepiped(
            height=size, width=size, thickness=size, nh=ncells, nw=ncells, nth=ncells)
    assert np.allclose(
            rp.faces_areas, [(size/(ncells-1))**2] * rp.nb_faces, rtol=1e-3)
    assert np.allclose(
            rp.faces_radiuses, [size/(ncells-1)*np.sqrt(2)/2] * rp.nb_faces, rtol=1e-3)
