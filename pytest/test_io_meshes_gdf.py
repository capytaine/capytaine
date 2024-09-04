"""Tests for the mesh import of mesh from WAMIT file format"""
from unittest import mock

import numpy as np
import pytest

import capytaine as cpt
from capytaine.io.mesh_loaders import load_GDF
from capytaine.io.mesh_writers import write_GDF


def parse_gdf(file_content):
    # Workaround to avoid actually writing/reading a file on disk
    with mock.patch('capytaine.io.mesh_loaders._check_file', lambda foo: None):
        with mock.patch('capytaine.io.mesh_loaders.open', mock.mock_open(read_data=file_content)):
            return load_GDF("mocked/filename")


def test_gdf_no_symmetry():
    mesh = parse_gdf("""Mesh description
        1.0  9.81
        0 0
        2
        0.0 0.0 -1.0
        1.0 0.0 -1.0
        1.0 1.0 -1.0
        0.0 1.0 -1.0
        1.0 0.0 -1.0
        2.0 0.0 -1.0
        2.0 1.0 -1.0
        1.0 1.0 -1.0
    """)
    assert isinstance(mesh, cpt.Mesh)
    assert mesh.nb_faces == 2


def test_gdf_no_symmetry_alternative_format():
    mesh = parse_gdf("""Mesh description
        1.0  9.81
        0 0
        2
        0.0 0.0 -1.0 1.0 0.0 -1.0 1.0 1.0 -1.0 0.0 1.0 -1.0
        1.0 0.0 -1.0 2.0 0.0 -1.0 2.0 1.0 -1.0 1.0 1.0 -1.0
    """)
    assert mesh.nb_faces == 2


def test_gdf_x_symmetry():
    mesh = parse_gdf("""Mesh description
        1.0  9.81
        1 0
        2
        0.0 0.0 -1.0
        1.0 0.0 -1.0
        1.0 1.0 -1.0
        0.0 1.0 -1.0
        1.0 0.0 -1.0
        2.0 0.0 -1.0
        2.0 1.0 -1.0
        1.0 1.0 -1.0
    """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert mesh.plane == cpt.yOz_Plane
    assert mesh.nb_faces == 4


def test_gdf_y_symmetry():
    mesh = parse_gdf("""Mesh description
        1.0  9.81
        0 1
        2
        0.0 0.0 -1.0
        1.0 0.0 -1.0
        1.0 1.0 -1.0
        0.0 1.0 -1.0
        1.0 0.0 -1.0
        2.0 0.0 -1.0
        2.0 1.0 -1.0
        1.0 1.0 -1.0
    """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert mesh.plane == cpt.xOz_Plane
    assert mesh.nb_faces == 4


def test_gdf_xy_symmetry():
    mesh = parse_gdf("""Mesh description
        1.0  9.81
        1 1
        2
        0.0 0.0 -1.0
        1.0 0.0 -1.0
        1.0 1.0 -1.0
        0.0 1.0 -1.0
        1.0 0.0 -1.0
        2.0 0.0 -1.0
        2.0 1.0 -1.0
        1.0 1.0 -1.0
    """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.ReflectionSymmetricMesh)
    assert mesh.nb_faces == 8


########################################################


def test_write_and_load_gdf(tmpdir):
    mesh_path = tmpdir.join("temp_mesh.gdf")

    original_mesh = cpt.mesh_horizontal_cylinder()
    write_GDF(str(mesh_path), original_mesh.vertices, original_mesh.faces, ulen=1, gravity=9.81, isx=0, isy=0)

    read_mesh = load_GDF(str(mesh_path))

    np.testing.assert_allclose(
        read_mesh.vertices[read_mesh.faces],
        original_mesh.vertices[original_mesh.faces],
        atol=1e-6
        )

    np.testing.assert_allclose(
        read_mesh.faces_areas,
        original_mesh.faces_areas
        )

    np.testing.assert_allclose(
        read_mesh.faces_normals,
        original_mesh.faces_normals,
        atol=1e-6
        )

    np.testing.assert_allclose(
        read_mesh.faces_centers,
        original_mesh.faces_centers
        )

    np.testing.assert_allclose(
        read_mesh.faces_radiuses,
        original_mesh.faces_radiuses
        )

    np.testing.assert_allclose(
        read_mesh.volume,
        original_mesh.volume
        )
