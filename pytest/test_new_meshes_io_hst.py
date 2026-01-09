# Copyright 2025 Mews Labs
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the mesh import of mesh from Hydrostar file format"""

from io import StringIO
import logging

import pytest
import numpy as np

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh


def parse_hst(content):
    return load_mesh(StringIO(content), file_format="hst")


def test_hst_type_0_triangle():
    """Load an HST mesh with PANEL TYPE 0 (no panel index) containing a quadrangle - verifies face center."""
    mesh = parse_hst(
        """
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """
    )
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1 / 2, 1 / 2, 0.0]]))


def test_hst_type_0_quadrangle():
    """Load an HST mesh with PANEL TYPE 0 containing a triangle (degenerate quad) - verifies face center."""
    mesh = parse_hst(
        """
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        2 3 4
        ENDPANEL
        """
    )
    assert mesh.vertices.shape == (3, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2 / 3, 2 / 3, 0.0]]))


def test_hst_type_1_quadrangle():
    """Load an HST mesh with PANEL TYPE 1 (with panel index) containing a quadrangle - verifies face center."""
    mesh = parse_hst(
        """
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 1 2 3 4
        ENDPANEL
        """
    )
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1 / 2, 1 / 2, 0.0]]))


def test_hst_type_1_triangle():
    """Load an HST mesh with PANEL TYPE 1 containing a triangle (degenerate quad) - verifies face center."""
    mesh = parse_hst(
        """
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 2 3 4
        ENDPANEL
        """
    )
    assert mesh.vertices.shape == (3, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2 / 3, 2 / 3, 0.0]]))


def test_hst_implicit_coordinate_numbering():
    """Load an HST mesh where vertex indices are implicit (not explicitly numbered) - verifies correct parsing."""
    mesh = parse_hst(
        """
        COORDINATES
        0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """
    )
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1 / 2, 1 / 2, 0.0]]))


def test_hst_coordinate_numbering_error():
    """Verify that HST parser raises ValueError when vertex numbering is incorrect or skips indices."""
    with pytest.raises(ValueError):
        parse_hst(
            """
            COORDINATES
            1 0.0 0.0 0.0
            2 1.0 0.0 0.0
            3 1.0 1.0 0.0
            5 0.0 1.0 0.0
            ENDCOORDINATES
            """
        )


def test_hst_symmetry_1():
    """Load an HST mesh with SYMMETRY 1 flag - expects a ReflectionSymmetricMesh with xOz plane."""
    mesh = parse_hst(
        """
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
            """
    )
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.plane == "xOz"
    assert mesh.nb_faces == 2  # Half mesh has 1 face, with symmetry = 2
    assert len(mesh.half.faces) == 1  # The actual half mesh


def test_hst_symmetry_2():
    """Load an HST mesh with SYMMETRY 2 flag - expects nested ReflectionSymmetricMesh with both planes."""
    mesh = parse_hst(
        """
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
            """
    )
    # Should return nested ReflectionSymmetricMesh
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert isinstance(mesh.half, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 4  # Quarter mesh has 1 face, with both symmetries = 4
    assert len(mesh.half.half.faces) == 1  # The actual quarter mesh


def test_hst_ignored_lines(caplog):
    """Verify that HST parser logs warnings for unrecognized keywords like USER and NBBODY."""
    with caplog.at_level(logging.WARNING):
        parse_hst(
            """
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
                """
        )
    assert "HST mesh reader ignored" in caplog.text
    assert "USER mancellin" in caplog.text
    assert "NBBODY 1" in caplog.text
