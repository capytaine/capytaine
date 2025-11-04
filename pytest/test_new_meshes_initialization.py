import numpy as np

import capytaine.new_meshes as meshes


def test_single_face_mesh():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2, 3]])
    mesh = meshes.Mesh(vertices=v, faces=f)
    assert mesh.nb_vertices == 4
    assert mesh.nb_faces == 1

    assert isinstance(mesh.faces_normals, np.ndarray)
    assert mesh.faces_normals.shape == (mesh.nb_faces, 3)
    assert np.allclose(mesh.faces_normals, np.array([[0.0, 0.0, 1.0]]))

    assert isinstance(mesh.faces_centers, np.ndarray)
    assert mesh.faces_centers.shape == (mesh.nb_faces, 3)
    assert np.allclose(mesh.faces_centers, np.array([0.5, 0.5, 0.0]))

    assert isinstance(mesh.faces_areas, np.ndarray)
    assert mesh.faces_areas.shape == (mesh.nb_faces,)
    assert np.allclose(mesh.faces_areas, np.array([1.0]))


def test_ignore_unused_vertices():
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 3.0, 0.0],
            [2.0, 3.0, 8.0],
        ]
    )
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 4


def test_ignore_duplicate_vertices():
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 4


def test_ignore_duplicate_faces():
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3], [1, 2, 3, 0]])
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 4


def test_normal_orientation():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    mesh1 = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    mesh2 = meshes.Mesh(vertices=v, faces=[[3, 2, 1, 0]])
    assert np.allclose(mesh1.faces_normals, -mesh2.faces_normals)


def test_areas_2():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    assert np.allclose(mesh.faces_areas, [2.0])


def test_areas_4():
    v = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    assert np.allclose(mesh.faces_areas, [4.0])


def test_centers():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 2.0, 0.0]])
    mesh = meshes.Mesh(vertices=v, faces=[[0, 1, 2, 3]])
    assert np.allclose(
        mesh.faces_centers, [0.44444, 0.77777, 0.0], atol=1e-3
    )  # Center of mass of the panel
    assert np.allclose(
        mesh.faces_vertices_centers, [0.5, 0.75, 0.0], atol=1e-3
    )  # Mean of the vertices


def test_mesh_with_triangle():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2, 2]])
    mesh = meshes.Mesh(vertices=v, faces=f)
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 3
    assert np.allclose(mesh.faces_normals, np.array([[0.0, 0.0, 1.0]]))


def test_from_list_of_faces():
    # Alternative way to define a mesh: array of shape (nb_faces, 4, 3)
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            [[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 1.0, 0.0], [2.0, 1.0, 0.0]],
        ]
    )
    assert mesh.nb_faces == 3
    assert mesh.nb_vertices == 8


def test_translation():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2, 3]])
    mesh1 = meshes.Mesh(vertices=v, faces=f)
    mesh2 = mesh1.translated(dx=1.0)
    assert mesh2.nb_faces == mesh1.nb_faces
    assert np.allclose(mesh2.faces_normals, mesh1.faces_normals)
    assert np.allclose(mesh2.faces_centers, np.array([1.5, 0.5, 0.0]))

    mesh3 = mesh1.translated(dz=1.0)
    assert mesh2.nb_faces == mesh1.nb_faces
    assert np.allclose(mesh2.faces_normals, mesh1.faces_normals)
    assert np.allclose(mesh3.faces_centers, np.array([0.5, 0.5, 1.0]))


def test_rotation():
    # Test rotations around x-axis
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2, 3]])
    mesh1 = meshes.Mesh(vertices=v, faces=f)

    mesh2 = mesh1.rotated(angle=np.pi / 2, axis="x")
    assert np.allclose(mesh2.faces_normals, np.array([[0.0, -1.0, 0.0]]))

    mesh3 = mesh1.rotated(angle=-np.pi / 2, axis="x")
    assert np.allclose(mesh3.faces_normals, np.array([[0.0, 1.0, 0.0]]))

    mesh4 = mesh1.rotated(angle=np.pi, axis="x")
    assert np.allclose(mesh4.faces_normals, np.array([[0.0, 0.0, -1.0]]))


def test_join_meshes():
    mesh1 = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ]
    )

    mesh2 = meshes.Mesh.from_list_of_faces(
        [
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ]
    )

    joined_mesh = mesh1 + mesh2
    assert joined_mesh.nb_faces == 2
    assert joined_mesh.nb_vertices == 6  # Duplicate vertices had been removed


def test_leading_count_column_uniform_quads():
    """Test that leading count column is properly stripped for uniform quad meshes."""
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
    )
    # Both faces are quads with leading count=4
    f = np.array([[4, 0, 1, 2, 3], [4, 1, 4, 5, 2]])
    mesh = meshes.Mesh(vertices=v, faces=f, auto_clean=False, auto_check=False)

    # The count column should be stripped
    assert mesh.faces.shape == (2, 4)
    assert mesh._faces == [[0, 1, 2, 3], [1, 4, 5, 2]]
    assert mesh.nb_faces == 2


def test_leading_count_column_mixed_triangles_quads():
    """Test that leading count column is properly stripped for mixed triangle/quad meshes."""
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
        ]
    )
    # First face: quad (4 vertices), Second face: triangle (3 vertices, padded)
    f = np.array([[4, 0, 1, 2, 3], [3, 1, 5, 4, 1]])
    mesh = meshes.Mesh(vertices=v, faces=f, auto_clean=False, auto_check=False)

    # The count column should be stripped
    assert mesh.faces.shape == (2, 4)
    assert mesh._faces == [[0, 1, 2, 3], [1, 5, 4, 1]]
    assert mesh.nb_faces == 2


def test_no_leading_count_column():
    """Test that faces without leading count column are processed correctly."""
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    # Normal face definition without count column
    f = np.array([[0, 1, 2, 3]])
    mesh = meshes.Mesh(vertices=v, faces=f, auto_clean=False, auto_check=False)

    # Should remain unchanged
    assert mesh.faces.shape == (1, 4)
    assert mesh._faces == [[0, 1, 2, 3]]
    assert mesh.nb_faces == 1


def test_leading_count_column_triangles_only():
    """Test that leading count column works for triangle-only meshes."""
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ]
    )
    # Two triangles with leading count=3 (padded to make uniform array)
    f = np.array([[3, 0, 1, 4, 4], [3, 1, 2, 4, 4]])
    mesh = meshes.Mesh(vertices=v, faces=f, auto_clean=False, auto_check=False)

    # The count column should be stripped
    assert mesh.faces.shape == (2, 4)
    assert mesh._faces == [[0, 1, 4, 4], [1, 2, 4, 4]]
    assert mesh.nb_faces == 2


def test_faces_as_list_no_count_column():
    """Test that list input is not affected by count column detection."""
    v = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    # List input with values that could be misinterpreted as count
    f = [[4, 0, 1, 2], [3, 1, 2, 3]]
    mesh = meshes.Mesh(vertices=v, faces=f, auto_clean=False, auto_check=False)

    # List input should pass through unchanged
    assert mesh._faces == [[4, 0, 1, 2], [3, 1, 2, 3]]
    assert mesh.nb_faces == 2
