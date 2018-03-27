#!/usr/bin/env python
# coding: utf-8

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import FloatingBody

from capytaine.meshes_collection import CollectionOfMeshes
from capytaine.symmetries import TranslationalSymmetry

from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.rectangle import Rectangle, OpenRectangularParallelepiped, RectangularParallelepiped


def test_collection_of_meshes():
    # Create some dummy meshes
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=np.float)
    dummy_meshes = [Mesh(vertices, [[0, 1, 2, 3]])]
    for i in range(3):
        dummy_meshes.append(dummy_meshes[0].copy())
        dummy_meshes[i+1].translate_z(i+1)

    # A first collection from a list
    coll = CollectionOfMeshes(dummy_meshes[:2])
    assert coll.nb_submeshes == 2
    assert coll.nb_vertices == 8
    assert coll.nb_faces == 2
    assert coll.submeshes[1].nb_faces == 1
    assert np.all(coll.faces_areas == np.asarray([1.0, 1.0]))
    assert np.all(coll.faces == np.asarray([[0, 1, 2, 3], [4, 5, 6, 7]]))
    assert np.all(coll.indices_of_mesh(1) == slice(1, 2))

    # A copy of the collection
    copy_coll = coll.copy()
    copy_coll.translate_x(1.0)  # Move
    assert copy_coll.nb_faces == 2
    assert np.all(copy_coll.vertices[:, 0] >= 1.0)  # Has moved

    # Another collection from an iterable
    other_coll = CollectionOfMeshes(iter(dummy_meshes))
    assert other_coll.nb_faces == 4
    assert np.all(other_coll.vertices[:, 0] <= 1.0)  # Did not move

    # A collection of collections
    big_coll = CollectionOfMeshes((copy_coll, other_coll))
    assert big_coll.nb_faces == 6
    assert big_coll.nb_submeshes == 2

    # Move one object in one of the sub-meshes
    copy_coll.submeshes[1].translate_x(1.0)
    assert big_coll.vertices[:, 0].max() == 3.0

    # Merging the big collection
    merged = big_coll.merge()
    assert isinstance(merged, Mesh)


def test_dof():
    nodes = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]])
    faces = np.array([[0, 1, 2, 3]])
    body = FloatingBody(Mesh(nodes, faces), name="one_face")
    assert body.dofs == {}

    body.add_translation_dof(direction=(1.0, 0.0, 0.0), name="1")
    assert body.dofs["1"] == np.array([0.0])

    body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    assert body.dofs["2"] == np.array([1.0])

    body.add_rotation_dof(axis_direction=(0.0, 0.0, 1.0), name="3")
    assert body.dofs["3"] == np.array([1.0])

    body.add_rotation_dof(axis_point=(0.5, 0, 0), axis_direction=(0.0, 0.0, 1.0), name="4")
    assert body.dofs["4"] == np.array([0.0])


def test_dof_name_inference():
    body = Sphere()
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    body.dofs['Surge_2'] = body.mesh.faces_normals @ (1, 0, 0)
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_1'])
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_2'])


def test_rectangle_generation():
    rec = Rectangle(size=(10, 10), resolution=(5, 2), clever=False, center=(0, 0, -5), name="test")
    assert rec.name == "test"
    assert isinstance(rec.mesh, Mesh)
    assert rec.mesh.name == "test_mesh"
    assert np.allclose(rec.size, np.array([10.0, 10.0]))
    assert np.allclose(rec.center, np.array([0, 0, -5.0]))
    assert rec.mesh.nb_faces == 10
    assert rec.mesh.nb_vertices == 18
    assert rec.area == 100

    assert np.all(rec.mesh.vertices[:, 0:2] <= 5.0)
    assert np.all(rec.mesh.vertices[:, 0:2] >= -5.0)
    assert np.all(rec.mesh.vertices[:, 2] <= 0.0)
    assert np.all(rec.mesh.vertices[:, 2] >= -10.0)

    # rec.show()

    rec2 = Rectangle(size=(10, 10), resolution=(5, 2), clever=True, center=(0, 0, -5), name="clever_test")
    assert isinstance(rec2.mesh, TranslationalSymmetry)
    assert rec2.mesh.nb_submeshes == 5
    assert rec2.mesh.nb_faces == 10

    merged_mesh = rec2.mesh.merge()
    assert merged_mesh.nb_vertices == 18
    assert merged_mesh.nb_faces == 10


def test_parallelepiped_generation():
    para = OpenRectangularParallelepiped(size=(5.0, 1.0, 3.0), center=(0, 0, -1.5),
                                         resolution=(5, 2, 3),
                                         clever=False, name="test")
    assert para.name == "test"
    assert isinstance(para.mesh, Mesh)
    assert para.mesh.name == "test_mesh"
    assert para.mesh.nb_faces == 42

    assert np.isclose(np.abs(para.mesh.vertices[:, 0]).max(), 2.5)
    assert np.isclose(np.abs(para.mesh.vertices[:, 1]).max(), 0.5)
    assert np.all(para.mesh.vertices[:, 2] <= 0.0)
    assert np.all(para.mesh.vertices[:, 2] >= -3.0)

    clever_para = OpenRectangularParallelepiped(size=(5.0, 1.0, 3.0), center=(0, 0, -1.5),
                                                resolution=(5, 2, 3),
                                                clever=True, name="clever_test")

    assert clever_para.mesh.nb_faces == 42
    assert isinstance(clever_para.mesh, CollectionOfMeshes)
    assert any(isinstance(mesh, TranslationalSymmetry)for mesh in clever_para.mesh.submeshes)

    full_para = RectangularParallelepiped(size=(5.0, 1.0, 3.0), center=(0, 0, -1.5),
                                          resolution=(5, 2, 3),
                                          clever=False, name="full_test")
    assert full_para.mesh.nb_faces == 62
    # full_para.show()

    clever_full_para = RectangularParallelepiped(size=(5.0, 1.0, 3.0), center=(0, 0, -1.5),
                                                 resolution=(5, 2, 3),
                                                 clever=True, name="clever_full_test")
    assert clever_full_para.mesh.nb_faces == 62
    # clever_full_para.show()
