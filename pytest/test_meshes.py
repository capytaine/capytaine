"""Tests for the mesh submodule: definition and transformation of base meshes."""

import pytest

import numpy as np
from numpy.linalg import norm
from numpy.typing import NDArray

import capytaine as cpt
from capytaine.meshes.properties import clustering

RNG = np.random.default_rng()


# Some meshes that will be used in the following tests.
test_mesh = cpt.Mesh(vertices=RNG.random((4, 3)), faces=[range(4)], name="test_mesh")
cylinder = cpt.mesh_horizontal_cylinder()
sphere = cpt.mesh_sphere(radius=1)


def test_mesh_initialization():
    """Test how the code checks the validity of the parameters."""
    with pytest.raises(AssertionError):
        cpt.Mesh(vertices=RNG.random((6, 3)), faces=[(0, 1, 2), (3, 4, 5)])

    with pytest.raises(AssertionError):
        cpt.Mesh(vertices=RNG.random((4, 3)), faces=[(0, 1, 2, -1)])

    with pytest.raises(AssertionError):
        cpt.Mesh(vertices=RNG.random((3, 3)), faces=[(0, 1, 2, 3)])


def test_vertices_and_faces_get_and_set():
    """Test get and set functions for vertices and faces."""
    # Get faces
    faces = cylinder.faces
    vertices = cylinder.vertices

    new_cylinder = cpt.Mesh(name="new_cylinder")

    # Set faces
    with pytest.raises(AssertionError):
        new_cylinder.faces = faces
        # You should set the vertices first.
    new_cylinder.vertices = vertices
    new_cylinder.faces = faces

    assert new_cylinder == cylinder


def test_mesh_naming():
    """Test how the mesh handle names and string representation."""
    # Test automatic naming
    dummy_mesh = cpt.Mesh()  # Automatically named something like mesh_1
    other_dummy_mesh = cpt.Mesh()  # Automatically named something like mesh_2
    assert dummy_mesh.name[:5] == "mesh_"
    assert other_dummy_mesh.name[:5] == "mesh_"
    assert int(dummy_mesh.name[5:]) + 1 == int(other_dummy_mesh.name[5:])


def test_as_set_of_faces():
    """Test the representation of the mesh as a set of faces.
    Also allows to define the equality of two meshes."""
    faces = cylinder.as_set_of_faces()
    assert isinstance(faces, frozenset)
    assert len(faces) == cylinder.nb_faces
    assert all(len(face) == 4 or len(face) == 3 for face in faces)  # Each face is represented by 3 or 4 points
    assert all(len(vertex) == 3 for face in faces for vertex in face)  # Each point is represented by 3 coordinates.

    assert cylinder == cylinder  # The equality is defined as the equality of the set of faces.
    assert cpt.Mesh.from_set_of_faces(faces) == cylinder  # The mesh can be reconstructed from a set of faces.


def mesh_from_set_of_faces():
    A, B, C, D = (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)

    # Two triangular faces
    faces = {frozenset({A, B, C}), frozenset({B, C, D})}
    mesh = cpt.Mesh.from_set_of_faces(faces)
    assert mesh.nb_vertices == 4
    assert mesh.nb_faces == 2
    assert mesh.is_triangle(0) and mesh.is_triangle(1)


def test_copy():
    """Test the copy and renaming of meshes."""
    assert test_mesh.copy() is not test_mesh
    assert test_mesh.copy() == test_mesh
    assert test_mesh.copy(name="copy_of_the_test_mesh").name == "copy_of_the_test_mesh"
    assert test_mesh.merged() is test_mesh


def test_shape_helper_functions():
    assert np.allclose(sphere.center_of_mass_of_nodes, (0.0, 0.0, 0.0))
    assert np.isclose(sphere.diameter_of_nodes, 2.0)
    assert cylinder.squared_axis_aligned_bbox == (-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)


def test_translate():
    assert cylinder.translated_x(0) == cylinder
    assert cylinder.translated_y(0) == cylinder
    assert cylinder.translated_z(0) == cylinder
    assert cylinder.translated([0, 0, 0]) == cylinder

    translated_sphere = sphere.translated([1, 2, 3])
    assert np.allclose(translated_sphere.center_of_mass_of_nodes, (1.0, 2.0, 3.0))
    assert np.isclose(translated_sphere.diameter_of_nodes, sphere.diameter_of_nodes)

def test_rotate():
    assert cylinder.rotated_x(0.0) == cylinder
    assert cylinder.rotated_y(0.0) == cylinder
    assert cylinder.rotated_z(0.0) == cylinder

    new_cylinder = cylinder.rotated_z(np.pi/2)
    assert np.allclose(new_cylinder.center_of_mass_of_nodes, cylinder.center_of_mass_of_nodes)
    assert new_cylinder.squared_axis_aligned_bbox == (-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)


def test_mirror():
    new_cylinder = cylinder.mirrored(cpt.Plane())
    cylinder.mirror(cpt.Plane())
    assert new_cylinder == cylinder


def test_symmetrized():
    from capytaine.meshes.symmetric import ReflectionSymmetricMesh
    sym = cylinder.merged().symmetrized(cpt.xOz_Plane)
    assert isinstance(sym, ReflectionSymmetricMesh)


def test_mesh_quality():
    # TODO: Could be improved...
    cylinder.triangulate_quadrangles()
    cylinder.merge_duplicates(atol=1e-5)
    cylinder.heal_mesh()

def test_heal_mesh_removes_degenerate_panels():
    vertices = np.array([(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)])
    faces = np.array([[0, 1, 2, 3], [1, 2, 2, 1]])
    mesh = cpt.Mesh(vertices, faces)
    assert mesh.nb_faces == 2
    mesh.heal_mesh()
    assert mesh.nb_faces == 1

def test_clipper():
    """Test clipping of mesh."""
    mesh = cpt.mesh_sphere(radius=5.0, resolution=(10, 5))
    aabb = mesh.axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, water_depth=np.inf)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, water_depth=1.0)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed

    # With CollectionOfMeshes (AxialSymmetry)
    mesh = cpt.mesh_sphere(radius=5.0, resolution=(10, 5), axial_symmetry=True)
    aabb = mesh.merged().axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, water_depth=np.inf)
    assert np.allclose(mesh.merged().axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, water_depth=1.0)
    assert np.allclose(mesh.merged().axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed

    # Check boundaries after clipping
    mesh = cpt.mesh_rectangle(size=(5,5), normal=(1,0,0))
    assert max([i[2] for i in mesh.immersed_part(free_surface=-1).vertices])<=-1
    assert max([i[2] for i in mesh.immersed_part(free_surface= 1).vertices])<= 1
    assert min([i[2] for i in mesh.immersed_part(free_surface=100, sea_bottom=-1).vertices])>=-1
    assert min([i[2] for i in mesh.immersed_part(free_surface=100, sea_bottom= 1).vertices])>= 1

    mesh = cpt.mesh_rectangle(size=(4,4), resolution=(1,1), normal=(1,0,0))
    tmp = list(mesh.clip(cpt.Plane(normal=(0,0.1,1),point=(0,0,-1)),inplace=False).vertices)
    tmp.sort(key=lambda x: x[2])
    tmp.sort(key=lambda x: x[1])
    assert np.allclose([i[2] for i in tmp], [-2, -0.8, -2, -1.2])


@pytest.mark.parametrize("size", [5, 6])
def test_clipper_indices(size):
    """Test clipped_mesh_faces_ids."""
    mesh = cpt.mesh_rectangle(size=(size, size), resolution=(size, size), center=(0, 0, 0))
    clipped_mesh = mesh.clipped(plane=cpt.Plane(point=(0, 0, 0), normal=(0, 0, 1)))
    faces_ids = clipped_mesh._clipping_data['faces_ids']

    assert clipped_mesh.nb_faces == len(faces_ids)
    assert all(norm(clipped_mesh.faces_centers[i] - mesh.faces_centers[face_id]) < 0.3
               for i, face_id in enumerate(faces_ids))


def test_clipper_corner_cases():
    mesh = sphere.translated_z(10.0)

    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, 1))
    clipped_mesh = mesh.clip(plane, inplace=False)
    assert clipped_mesh == cpt.Mesh(None, None)  # Empty mesh

    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, -1))
    clipped_mesh = mesh.clip(plane, inplace=False)
    assert clipped_mesh == mesh  # Unchanged mesh

    # Two distinct bodies
    two_spheres = cpt.Mesh.join_meshes(sphere.translated_z(10.0), sphere.translated_z(-10.0))
    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, -1))
    one_sphere_remaining = two_spheres.clip(plane, inplace=False)
    assert one_sphere_remaining == sphere.translated_z(10.0)


def test_clipper_tolerance():
    mesh = cpt.mesh_vertical_cylinder(length=10.001, center=(0, 0, -5))
    mesh = mesh.immersed_part()
    np.testing.assert_allclose(mesh.vertices[:, 2].max(), 0.0, atol=1e-12)


def test_extract_one_face():
    i = 2
    one_face = sphere.extract_one_face(i)
    assert np.all(one_face.faces_centers[0] == sphere.faces_centers[i])


##########################
#  Connected components  #
##########################

def test_connected_components_of_empty_mesh():
    from capytaine.meshes.properties import connected_components
    mesh = cpt.Mesh()
    cc = connected_components(mesh)
    assert len(cc) == 0


def test_connected_components_of_sphere():
    from capytaine.meshes.properties import connected_components
    mesh = cpt.mesh_sphere()
    cc = connected_components(mesh)
    assert len(cc) == 1


def test_connected_components_of_two_spheres():
    from capytaine.meshes.properties import connected_components
    mesh = cpt.mesh_sphere() + cpt.mesh_sphere(center=(5, 0, 0))
    cc = connected_components(mesh)
    assert len(cc) == 2


def test_connected_components_of_open_sphere():
    from capytaine.meshes.properties import connected_components
    mesh = cpt.mesh_sphere().immersed_part()
    cc = connected_components(mesh)
    assert len(cc) == 1


def test_connected_components_of_torus():
    from capytaine.io.meshio import load_from_meshio
    from capytaine.meshes.properties import connected_components
    pygmsh = pytest.importorskip("pygmsh")
    with pygmsh.occ.Geometry() as geom:
        geom.add_torus((0, 0, 0), 2.0, 0.5, mesh_size=0.3)
        gmsh_mesh = geom.generate_mesh(dim=2)
    mesh = load_from_meshio(gmsh_mesh).heal_mesh()
    cc = connected_components(mesh)
    assert len(cc) == 1


######################################
#  Connected components at waterline #
######################################

def test_connected_components_at_waterline_of_sphere():
    from capytaine.meshes.properties import connected_components_of_waterline
    mesh = cpt.mesh_sphere()
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 1


def test_connected_components_at_waterline_of_two_spheres():
    from capytaine.meshes.properties import connected_components_of_waterline
    mesh = cpt.mesh_sphere() + cpt.mesh_sphere(center=(5, 0, 0))
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 2


def test_connected_components_at_waterline_of_open_sphere():
    from capytaine.meshes.properties import connected_components_of_waterline
    mesh = cpt.mesh_sphere().immersed_part()
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 1


def test_connected_components_at_waterline_of_immersed_sphere():
    from capytaine.meshes.properties import connected_components_of_waterline
    mesh = cpt.mesh_sphere(center=(0, 0, -5))
    cc = connected_components_of_waterline(mesh, z=0.0)
    assert len(cc) == 0


def test_connected_components_at_waterline_of_torus():
    from capytaine.io.meshio import load_from_meshio
    from capytaine.meshes.properties import connected_components_of_waterline
    pygmsh = pytest.importorskip("pygmsh")
    with pygmsh.occ.Geometry() as geom:
        geom.add_torus((0, 0, 0), 2.0, 0.5, mesh_size=0.3)
        gmsh_mesh = geom.generate_mesh(dim=2)
    mesh = load_from_meshio(gmsh_mesh).heal_mesh()
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 2

####################
#  Lid extraction  #
####################

def test_extract_lid():
    mesh = cpt.mesh_vertical_cylinder(center=(0, 0, -2), length=4, resolution=(2, 8, 10))
    hull_mesh, lid_mesh = mesh.extract_lid()
    assert lid_mesh.nb_faces == 2*8


###################################
#  Connected vertices clustering  #
###################################

@pytest.mark.parametrize("faces", [
    np.array([[1, 2], [1, 3], [2, 3], [4, 5], [5, 6], [7, 8], [3, 9]]),
    np.array([[0, 1, 2, 3], [1, 2, 6, 7]]),
])
def test_vertice_clustering(faces: NDArray[np.integer]):
    """Test the clustering algorithm for connected faces & vertices."""
    # Legacy way to cluster
    vertices_components: set[frozenset[int]] = set()
    for set_of_v_in_face in map(frozenset, faces):
        intersecting_components = [c for c in vertices_components if len(c.intersection(set_of_v_in_face)) > 0]
        if len(intersecting_components) == 0:
            vertices_components.add(set_of_v_in_face)
        else:
            for c in intersecting_components:
                vertices_components.remove(c)
            vertices_components.add(frozenset.union(set_of_v_in_face, *intersecting_components))
    # Vectorized clustering
    vert_groups = clustering(faces)
    assert {frozenset({*group}) for group in vert_groups} == vertices_components
