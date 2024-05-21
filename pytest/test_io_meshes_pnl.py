"""Tests for the mesh import of mesh from HAMS file format"""
from unittest import mock

import numpy as np

import capytaine as cpt
from capytaine.io.mesh_loaders import load_PNL
from capytaine.io.mesh_writers import write_PNL


def parse_pnl(file_content):
    # Workaround to avoid actually writing/reading a file on disk
    with mock.patch('capytaine.io.mesh_loaders._check_file', lambda foo: None):
        with mock.patch('capytaine.io.mesh_loaders.open', mock.mock_open(read_data=file_content)):
            return load_PNL("mocked/filename")


def test_pnl_no_symmetry():
    mesh = parse_pnl(
            """ --------------Hull Mesh File---------------

            # Number of Panels, Nodes, X-Symmetry and Y-Symmetry
            2        8           0           0

            #Start Definition of Node Coordinates     ! node_number   x   y   z
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            #End Definition of Node Coordinates

            #Start Definition of Node Relations
            1 4  1 2 3 4
            2 4  5 6 7 8
            #End Definition of Node Relations
            """)
    assert isinstance(mesh, cpt.Mesh)
    assert mesh.nb_faces == 2


def test_pnl_triangle():
    mesh = parse_pnl(
            """ --------------Hull Mesh File---------------

            # Number of Panels, Nodes, X-Symmetry and Y-Symmetry
            2        7           0           0

            #Start Definition of Node Coordinates     ! node_number   x   y   z
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            #End Definition of Node Coordinates

            #Start Definition of Node Relations
            1 4  1 2 3 4
            2 3  5 6 7
            #End Definition of Node Relations
            """)
    assert isinstance(mesh, cpt.Mesh)
    assert mesh.nb_faces == 2


def test_pnl_x_symmetry():
    mesh = parse_pnl(
            """ --------------Hull Mesh File---------------

            # Number of Panels, Nodes, X-Symmetry and Y-Symmetry
            2        8           1           0

            #Start Definition of Node Coordinates     ! node_number   x   y   z
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            #End Definition of Node Coordinates

            #Start Definition of Node Relations
            1 4  1 2 3 4
            2 4  5 6 7 8
            #End Definition of Node Relations
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert mesh.plane == cpt.yOz_Plane
    assert mesh.nb_faces == 4


def test_pnl_y_symmetry():
    mesh = parse_pnl(
            """ --------------Hull Mesh File---------------

            # Number of Panels, Nodes, X-Symmetry and Y-Symmetry
            2        8           0           1

            #Start Definition of Node Coordinates     ! node_number   x   y   z
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            #End Definition of Node Coordinates

            #Start Definition of Node Relations
            1 4  1 2 3 4
            2 4  5 6 7 8
            #End Definition of Node Relations
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert mesh.plane == cpt.xOz_Plane
    assert mesh.nb_faces == 4


def test_pnl_xy_symmetry():
    mesh = parse_pnl(
            """ --------------Hull Mesh File---------------

            # Number of Panels, Nodes, X-Symmetry and Y-Symmetry
            2        8           1           1

            #Start Definition of Node Coordinates     ! node_number   x   y   z
            1 0.0 0.0 -1.0
            2 1.0 0.0 -1.0
            3 1.0 1.0 -1.0
            4 0.0 1.0 -1.0
            5 1.0 0.0 -1.0
            6 2.0 0.0 -1.0
            7 2.0 1.0 -1.0
            8 1.0 1.0 -1.0
            #End Definition of Node Coordinates

            #Start Definition of Node Relations
            1 4  1 2 3 4
            2 4  5 6 7 8
            #End Definition of Node Relations
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.ReflectionSymmetricMesh)
    assert mesh.nb_faces == 8


###################################################


def test_pnl_roundtrip(tmpdir):
    mesh_path = tmpdir.join("temp_mesh.pnl")
    mesh = cpt.mesh_sphere()
    write_PNL(mesh_path, mesh.vertices, mesh.faces)
    reloaded_mesh = cpt.load_mesh(mesh_path)
    np.testing.assert_equal(mesh.faces, reloaded_mesh.faces)
    np.testing.assert_allclose(mesh.vertices, reloaded_mesh.vertices, atol=1e-5)
