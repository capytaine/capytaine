"""Tests for the mesh import of mesh from Hydrostar file format"""
from unittest import mock
import logging

import pytest
import numpy as np

import capytaine as cpt
from capytaine.io.mesh_loaders import load_HST


def parse_hst(file_content):
    # Workaround to avoid actually writing/reading a file on disk
    with mock.patch('capytaine.io.mesh_loaders._check_file', lambda foo: None):
        with mock.patch('capytaine.io.mesh_loaders.open', mock.mock_open(read_data=file_content)):
            return load_HST("mocked/filename")


def test_hst_type_0_triangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_type_0_quadrangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2/3, 2/3, 0.0]]))


def test_hst_type_1_quadrangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_type_1_triangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2/3, 2/3, 0.0]]))


def test_hst_implicit_coordinate_numbering():
    mesh = parse_hst("""
        COORDINATES
        0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_coordinate_numbering_error():
    with pytest.raises(ValueError):
        parse_hst("""
            COORDINATES
            1 0.0 0.0 0.0
            2 1.0 0.0 0.0
            3 1.0 1.0 0.0
            5 0.0 1.0 0.0
            ENDCOORDINATES
            """)

def test_hst_symmetry_1():
    mesh = parse_hst("""
        SYMMETRY 1
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.Mesh)
    assert mesh.merged().nb_faces == 2


def test_hst_symmetry_2():
    mesh = parse_hst("""
        SYMMETRY 2
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half.half, cpt.Mesh)
    assert mesh.merged().nb_faces == 4


def test_hst_ignored_lines(caplog):
    with caplog.at_level(logging.WARNING):
        parse_hst("""
            USER mancellin
            NBBODY 1
            COORDINATES
            1 0.0 0.0 0.0
            2 1.0 0.0 0.0
            3 1.0 1.0 0.0
            4 0.0 1.0 0.0
            ENDCOORDINATES
            PANEL TYPE 0
            1 2 3 4
            ENDPANEL
                """)
    assert "HST mesh reader ignored" in caplog.text
    assert "USER mancellin" in caplog.text
    assert "NBBODY 1" in caplog.text
