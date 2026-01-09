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

"""Tests for the mesh import of mesh from HAMS file format"""

from io import StringIO

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh


def parse_pnl(content):
    return load_mesh(StringIO(content), file_format="pnl")


def test_pnl_no_symmetry():
    """Load a PNL mesh without symmetry flags - expects a basic Mesh with 2 quadrangular faces."""
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
            """
    )
    assert mesh.nb_faces == 2
    assert isinstance(mesh, Mesh) and not isinstance(mesh, ReflectionSymmetricMesh)


def test_pnl_triangle():
    """Load a PNL mesh containing both quadrangles and triangles - expects 2 faces total."""
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
            """
    )
    assert mesh.nb_faces == 2


def test_pnl_x_symmetry():
    """Load a PNL mesh with X-symmetry flag - expects ReflectionSymmetricMesh with xOz plane."""
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
            """
    )
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.plane == "yOz"
    assert mesh.nb_faces == 4  # Half mesh has 2 faces, with symmetry = 4
    assert len(mesh.half.faces) == 2


def test_pnl_y_symmetry():
    """Load a PNL mesh with Y-symmetry flag - expects ReflectionSymmetricMesh with yOz plane."""
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
            """
    )
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.plane == "xOz"
    assert mesh.nb_faces == 4
    assert len(mesh.half.faces) == 2


def test_pnl_xy_symmetry():
    """Load a PNL mesh with both X and Y symmetry flags - expects nested ReflectionSymmetricMesh."""
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
            """
    )
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert isinstance(mesh.half, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 8  # Quarter mesh has 2 faces, with both symmetries = 8
    assert len(mesh.half.half.faces) == 2
