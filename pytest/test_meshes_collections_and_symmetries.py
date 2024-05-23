"""Tests for the mesh submodule: definition and transformation of collections of meshes and symmetric meshes."""

import pytest
import numpy as np
import capytaine as cpt


def test_collection_of_meshes():
    # Create some dummy meshes
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    dummy_meshes = [cpt.Mesh(vertices, [[0, 1, 2, 3]])]
    for i in range(3):
        dummy_meshes.append(dummy_meshes[0].copy())
        dummy_meshes[i+1].translate_z(i+1)

    # A first collection from a list
    coll = cpt.CollectionOfMeshes(dummy_meshes[:2])
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
    other_coll = cpt.CollectionOfMeshes(iter(dummy_meshes))
    assert other_coll.nb_faces == 4
    assert np.all(other_coll.vertices[:, 0] <= 1.0)  # Did not move

    # A collection of collections
    big_coll = cpt.CollectionOfMeshes((copy_coll, other_coll))
    assert big_coll.nb_faces == 6
    assert big_coll.nb_submeshes == 2

    # Move one object in one of the sub-meshes
    copy_coll[1].translate_x(1.0)
    assert big_coll.vertices[:, 0].max() == 3.0

    # Merging the big collection
    merged = big_coll.merged()
    assert isinstance(merged, cpt.Mesh)


def test_collection():
    sphere = cpt.mesh_sphere(name="foo", center=(0, 0, -2))
    other_sphere = cpt.mesh_sphere(name="bar", center=(0, 0, 2))

    coll = cpt.CollectionOfMeshes([sphere, other_sphere], name="baz")
    coll2 = cpt.CollectionOfMeshes([sphere, other_sphere])

    assert coll == coll2
    assert hash(coll) == hash(coll2)

    assert coll[0] == coll2[0]

    coll.tree_view()

    clipped_coll = coll.keep_immersed_part(inplace=False)
    assert len(clipped_coll) == 1
    # assert clipped_coll.merged() == sphere.merged()

    assert np.allclose(sphere.center_of_mass_of_nodes, sphere.merged().center_of_mass_of_nodes)
    assert np.allclose(sphere.diameter_of_nodes, sphere.merged().diameter_of_nodes)


def test_join_reflection_symmetric_meshes():
    cylinder_1 = cpt.mesh_horizontal_cylinder(center=(0, 0, -2), reflection_symmetry=True)
    cylinder_2 = cpt.mesh_horizontal_cylinder(center=(0, 0, -4), reflection_symmetry=True)
    assert cylinder_1.plane == cylinder_2.plane
    both = cylinder_1.join_meshes(cylinder_2)
    assert isinstance(both, cpt.ReflectionSymmetricMesh)


def test_join_axisymmetric_disks():
    """Given two axisymmetric meshes with the same axis, build a new axisymmetric mesh combining the two."""
    disk1 = cpt.mesh_disk(radius=1.0, center=(-1, 0, 0), resolution=(6, 6), normal=(1, 0, 0), axial_symmetry=True)
    disk2 = cpt.mesh_disk(radius=2.0, center=(1, 0, 0), resolution=(8, 6), normal=(1, 0, 0), axial_symmetry=True)
    joined = disk1.join_meshes(disk2, name="two_disks")
    assert isinstance(joined, cpt.AxialSymmetricMesh)
    joined.tree_view()

    disk3 = cpt.mesh_disk(radius=1.0, center=(0, 0, 0), resolution=(6, 4), axial_symmetry=True)
    with pytest.raises(AssertionError):
        disk1.join_meshes(disk3)


def test_join_translational_cylinders():
    """Given two meshes with the same translation symmetry, join them into a single mesh with the same symmetry."""
    params = dict(length=10.0, reflection_symmetry=False, translation_symmetry=True, resolution=(0, 10, 10))
    mesh1 = cpt.mesh_horizontal_cylinder(radius=1.0, center=(0, 5, -5), **params)
    mesh2 = cpt.mesh_horizontal_cylinder(radius=2.0, center=(0, -5, -5), **params)
    joined = mesh1.join_meshes(mesh2)
    assert isinstance(joined, cpt.TranslationalSymmetricMesh)


def test_mesh_splitting():
    mesh = cpt.mesh_sphere()

    splitted_mesh = mesh.sliced_by_plane(cpt.xOz_Plane.translated_y(0.5))
    assert isinstance(splitted_mesh, cpt.CollectionOfMeshes)
    assert splitted_mesh.merged() == mesh

    twice_splitted_mesh = splitted_mesh.sliced_by_plane(cpt.yOz_Plane)
    assert isinstance(twice_splitted_mesh[0], cpt.CollectionOfMeshes)
    assert isinstance(twice_splitted_mesh[1], cpt.CollectionOfMeshes)
    assert twice_splitted_mesh.merged() == mesh


def test_extract_one_face():
    sphere = cpt.mesh_sphere(axial_symmetry=True)
    assert sphere.submesh_containing_face(0) == (0, 0)
    assert sphere.submesh_containing_face(1) == (0, 1)
    assert sphere.submesh_containing_face(sphere[0].nb_faces-1) == (0, sphere[0].nb_faces-1)
    assert sphere.submesh_containing_face(sphere[0].nb_faces) == (1, 0)
    assert sphere.submesh_containing_face(sphere[0].nb_faces+1) == (1, 1)

    for i in [5, sphere[0].nb_faces+5]:
        assert np.allclose(sphere.extract_one_face(i).faces_centers[0], sphere.faces_centers[i])


def test_immersed_part():
    mesh = cpt.mesh_horizontal_cylinder(reflection_symmetry=True)
    assert mesh.immersed_part().merged() == mesh.merged().immersed_part()


def test_path_to_leaf():
    sphere_1 = cpt.mesh_sphere()
    sphere_2 = cpt.mesh_sphere()
    sphere_3 = cpt.mesh_sphere()
    sphere_4 = cpt.mesh_sphere()
    coll = cpt.CollectionOfMeshes([
        cpt.CollectionOfMeshes([sphere_1, sphere_2]),
        sphere_3, sphere_4])
    assert sphere_1.path_to_leaf() == [[]]
    assert coll.path_to_leaf() == [[0, 0], [0, 1], [1], [2]]
