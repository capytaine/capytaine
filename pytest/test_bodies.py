#!/usr/bin/env python
# coding: utf-8

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import FloatingBody
from capytaine.meshes_collection import CollectionOfMeshes


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


# def test_collection():
    # body_1 = generate_sphere()
    # body_1.name = 'body_1'
    # body_1.dofs['Dummy'] = body_1.mesh.faces_normals @ (0, 0, 0)
    # body_1.dofs['Heave'] = body_1.mesh.faces_normals @ (0, 0, 1)
    # assert isinstance(body_1, FloatingBody)
    # assert isinstance(body_1.as_FloatingBody(), FloatingBody)
    #
    # body_2 = generate_sphere(z0=-5.0)
    # body_2.name = 'body_2'
    # body_2.dofs['Surge'] = body_2.mesh.faces_normals @ (1, 0, 0)
    # coll = body_1 + body_2
    # assert isinstance(coll, CollectionOfFloatingBodies)
    # assert isinstance(coll.as_FloatingBody(), FloatingBody)
    # assert coll.name == 'CollectionOfFloatingBodies([body_1, body_2])'
    #
    # assert np.all(coll.mesh.faces_areas[coll.indices_of_body(0)] == body_1.mesh.faces_areas)
    # assert np.all(coll.mesh.faces_areas[coll.indices_of_body(1)] == body_2.mesh.faces_areas)
    #
    # # Test dofs
    # assert len(coll.dofs) == 3
    # for name, dof in coll.dofs.items():
    #     assert name in ['body_1_Dummy', 'body_1_Heave', 'body_2_Surge']
    #     assert dof.shape[0] == coll.mesh.faces_normals.shape[0]
    #     assert np.all(dof[:dof.shape[0]//2] == np.zeros(coll.mesh.faces_normals.shape[0]//2)) \
    #             or np.all(dof[dof.shape[0]//2:] == np.zeros(coll.mesh.faces_normals.shape[0]//2))
    #
    # body_3 = generate_sphere(z0=-10.0)
    # body_3.name = 'body_3'
    # coll = body_1 + body_2 + body_3
    # assert isinstance(coll, CollectionOfFloatingBodies)
    # assert coll.name == 'CollectionOfFloatingBodies([body_1, body_2, body_3])'
    # assert body_1 in coll.subbodies
    # assert body_2 in coll.subbodies
    # assert body_3 in coll.subbodies
    #
    # assert coll.mesh.nb_vertices == coll.as_FloatingBody().mesh.nb_vertices == generate_sphere().mesh.nb_vertices*3
    # assert coll.mesh.nb_faces == coll.as_FloatingBody().mesh.nb_faces == generate_sphere().mesh.nb_faces*3
    #
    # # Test the merging of identical vertices
    # assert (generate_sphere() + generate_sphere()).as_FloatingBody().mesh.nb_vertices == generate_sphere().mesh.nb_vertices

#
# def test_symmetric_bodies():
#     half_sphere = generate_half_sphere(nphi=5)
#     half_sphere.name = 'half_sphere'
#     full_sphere = ReflectionSymmetry(half_sphere, xOz_Plane)
#     assert isinstance(full_sphere, CollectionOfFloatingBodies)
#     assert half_sphere in full_sphere.subbodies
#     assert full_sphere.as_FloatingBody().mesh.nb_vertices == generate_sphere(nphi=10).mesh.nb_vertices
#
#     other_sphere = generate_sphere(z0=-5.0)
#     coll = full_sphere + other_sphere
#     assert full_sphere in coll.subbodies

