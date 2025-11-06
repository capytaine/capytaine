import logging

import pytest
import numpy as np

import capytaine.new_meshes as meshes
from capytaine.new_meshes.quality import (
    indices_of_non_convex_faces,
    indices_of_non_coplanar_faces,
)

def test_nan_in_vertices():
    with pytest.raises(ValueError):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [[0.0, 0.0, 0.0], [1.0, np.nan, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        )

def test_invalid_indices():
    with pytest.raises(ValueError):
        v = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        meshes.Mesh(vertices=v, faces=[[0, 1, 2, 8]])

def test_invalid_indices():
    with pytest.raises(ValueError):
        v = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        meshes.Mesh(vertices=v, faces=[[0, 1, -1, 2]])

def test_coplanar_quads(caplog):
    """All quadrilaterals are coplanar → no indices returned."""
    with caplog.at_level(logging.WARNING):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        )
    assert caplog.text == ""
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == []


def test_one_non_coplanar_quad(caplog):
    """First quadrilateral is non-coplanar → index [0]."""
    with caplog.at_level(logging.WARNING):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        )
    assert "non-coplanar" in caplog.text
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == [0]


def test_two_non_coplanar_quads(caplog):
    """Both quadrilaterals are non-coplanar → indices [0, 1]."""
    with caplog.at_level(logging.WARNING):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.1], [1.0, 1.0, 0.0]],
            ]
        )
    assert "non-coplanar" in caplog.text
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == [0, 1]


def test_convex_quads_are_valid(caplog):
    """Convex quadrilaterals should not be flagged as non-convex."""
    with caplog.at_level(logging.WARNING):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        )
    assert caplog.text == ""
    assert indices_of_non_convex_faces(mesh.vertices, mesh.faces) == []


def test_non_convex_quad_is_detected(caplog):
    """A quadrilateral with crossed diagonals (vertex order issue) should be flagged."""
    with caplog.at_level(logging.WARNING):
        mesh = meshes.Mesh.from_list_of_faces(
            [
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],  # swapped vertices → concave / self-intersecting
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
            ]
        )
    assert "non-convex" in caplog.text
    assert indices_of_non_convex_faces(mesh.vertices, mesh.faces) == [0]
