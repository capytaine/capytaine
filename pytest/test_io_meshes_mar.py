"""Tests for the mesh import of mesh from NEMOH file format"""
from unittest import mock

import capytaine as cpt
from capytaine.io.mesh_loaders import load_MAR


def parse_mar(file_content):
    # Workaround to avoid actually writing/reading a file on disk
    with mock.patch('capytaine.io.mesh_loaders._check_file', lambda foo: None):
        with mock.patch('capytaine.io.mesh_loaders.open', mock.mock_open(read_data=file_content)):
            return load_MAR("mocked/filename")


def test_mar_no_symmetry():
    mesh = parse_mar(
            """2 0
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            0 0 0 0
            1 2 3 4
            5 6 7 8
            0 0 0 0
            """)
    assert isinstance(mesh, cpt.Mesh)
    assert mesh.nb_faces == 2


def test_mar_symmetry():
    mesh = parse_mar(
            """2 1
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            0 0 0 0
            1 2 3 4
            5 6 7 8
            0 0 0 0
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert mesh.plane == cpt.xOz_Plane
    assert mesh.nb_faces == 4
