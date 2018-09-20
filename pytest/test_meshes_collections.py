import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.mesh.meshes_collection import CollectionOfMeshes
from capytaine.geometric_bodies import Sphere


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
    assert coll[1].nb_faces == 1
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
    copy_coll[1].translate_x(1.0)
    assert big_coll.vertices[:, 0].max() == 3.0

    # Merging the big collection
    merged = big_coll.merge()
    assert isinstance(merged, Mesh)


def test_collection():
    sphere = Sphere(name="foo", center=(0, 0, -2)).mesh
    other_sphere = Sphere(name="bar", center=(0, 0, 2)).mesh

    coll = CollectionOfMeshes([sphere, other_sphere], name="baz")
    assert str(coll) == "baz"
    assert repr(coll) == "CollectionOfMeshes(('foo_mesh', 'bar_mesh'), name=baz)"

    coll2 = CollectionOfMeshes([sphere, other_sphere])
    assert repr(coll2) == "CollectionOfMeshes('foo_mesh', 'bar_mesh')"
    assert str(coll2) == repr(coll2)

    assert coll == coll2
    assert hash(coll) == hash(coll2)

    assert coll[0] == coll2[0]

    coll.tree_view()

    clipped_coll = coll.keep_immersed_part(inplace=False)
    assert len(clipped_coll) == 1
    merged = clipped_coll.merge()
    assert merged == sphere.merge()
