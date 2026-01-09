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

"""Tests for the mesh import of mesh from AQUADYN and NEMOH's file format"""

import os
import gzip
from io import StringIO

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh


def parse_mar(content):
    return load_mesh(StringIO(content), file_format="mar")


def test_mar_no_symmetry():
    """Load a MAR mesh without symmetry flags - expects a basic Mesh with 2 faces."""
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
            """
    )
    assert mesh.nb_faces == 2
    assert isinstance(mesh, Mesh) and not isinstance(mesh, ReflectionSymmetricMesh)


def test_mar_symmetry():
    """Load a MAR mesh with xOz symmetry flag - should be a ReflectionSymmetricMesh."""
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
            """
    )
    # With xOz symmetry, the mesh is reflected, so total faces = 2 * 2 = 4
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 4
    assert mesh.plane == "xOz"


def test_symmetric_mar_from_path():
    mesh = load_mesh(os.path.join(
        os.path.dirname(__file__),
        "Nemoh_verification_cases/Cylinder/Cylinder.dat"
    ), file_format="nemoh")
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 2*300


def test_symmetric_mar_from_file():
    with open(os.path.join(
        os.path.dirname(__file__),
        "Nemoh_verification_cases/Cylinder/Cylinder.dat"
        ), 'r') as f:
        mesh = load_mesh(f, file_format="nemoh")
    assert isinstance(mesh, ReflectionSymmetricMesh)
    assert mesh.nb_faces == 2*300


def test_mar_from_path():
    mesh = load_mesh(os.path.join(
        os.path.dirname(__file__),
        "Nemoh_verification_cases/NonSymmetrical/NonSymmetrical.dat"
    ), file_format="nemoh")
    assert isinstance(mesh, Mesh)
    assert mesh.nb_faces == 280


def test_mar_from_file():
    with open(os.path.join(
        os.path.dirname(__file__),
        "Nemoh_verification_cases/NonSymmetrical/NonSymmetrical.dat"
        ), 'r') as f:
        mesh = load_mesh(f, file_format="nemoh")
    assert isinstance(mesh, Mesh)
    assert mesh.nb_faces == 280


def test_mar_from_compressed_file():
    path = os.path.join(
        os.path.dirname(__file__),
        "mesh_files_examples/boat_200.mar.gz"
        )
    with gzip.open(path, 'rt') as f:
        mesh = load_mesh(f, file_format="nemoh")
    assert mesh.nb_faces == 500
