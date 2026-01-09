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

"""Tests for the mesh import of mesh from WAMIT file format"""

from io import StringIO

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh


def parse_gdf(content):
    return load_mesh(StringIO(content), file_format="gdf")


def test_gdf_no_symmetry():
    """Load a GDF mesh without symmetry flags - expects a basic Mesh with 2 faces."""
    mesh = parse_gdf(
        """Mesh description
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
    """
    )
    assert mesh.nb_faces == 2
    assert isinstance(mesh, Mesh)
    assert not isinstance(mesh, ReflectionSymmetricMesh)


def test_gdf_no_symmetry_alternative_format():
    """Load a GDF mesh with all face coordinates on a single line - alternative formatting."""
    mesh = parse_gdf(
        """Mesh description
        1.0  9.81
        0 0
        2
        0.0 0.0 -1.0 1.0 0.0 -1.0 1.0 1.0 -1.0 0.0 1.0 -1.0
        1.0 0.0 -1.0 2.0 0.0 -1.0 2.0 1.0 -1.0 1.0 1.0 -1.0
    """
    )
    assert mesh.nb_faces == 2


def test_gdf_x_symmetry():
    """Load a GDF mesh with X-symmetry flag (isx=1) - should be a ReflectionSymmetricMesh with yOz plane."""

    mesh = parse_gdf(
        """Mesh description
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
    """
    )
    # With yOz symmetry (x-symmetry), total faces = 2 * 2 = 4
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 4
    assert mesh.plane == "yOz"


def test_gdf_y_symmetry():
    """Load a GDF mesh with Y-symmetry flag (isy=1) - should be a ReflectionSymmetricMesh with xOz plane."""

    mesh = parse_gdf(
        """Mesh description
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
    """
    )
    # With xOz symmetry (y-symmetry), total faces = 2 * 2 = 4
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 4
    assert mesh.plane == "xOz"


def test_gdf_xy_symmetry():
    """Load a GDF mesh with both X and Y symmetry flags - should be a nested ReflectionSymmetricMesh."""

    mesh = parse_gdf(
        """Mesh description
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
    """
    )
    # With both symmetries (quarter mesh), total faces = 2 * 2 * 2 = 8
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert isinstance(mesh.half, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 8
    assert mesh.plane == "yOz"  # Outer symmetry
    assert mesh.half.plane == "xOz"  # Inner symmetry
