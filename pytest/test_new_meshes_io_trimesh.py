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

import os
import lzma
import pytest

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh

def test_load_directly_from_trimesh(tmp_path):
    """
    Test loading a simple triangle mesh from an .obj file using trimesh.
    """
    pytest.importorskip("trimesh", reason="trimesh not installed, test skipped")

    import trimesh
    trimesh_mesh = trimesh.Trimesh(
        vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], faces=[[0, 1, 2]]
    )
    mesh = load_mesh(trimesh_mesh)
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 3


def test_load_from_trimesh(tmp_path):
    """
    Test loading a simple triangle mesh from an .obj file using trimesh.
    """
    pytest.importorskip("trimesh", reason="trimesh not installed, test skipped")

    import trimesh
    trimesh_mesh = trimesh.Trimesh(
        vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], faces=[[0, 1, 2]]
    )

    test_file = tmp_path / "dummy.glb"
    trimesh_mesh.export(test_file)

    mesh = load_mesh(test_file, "glb")
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 3


def test_export_to_trimesh():
    pytest.importorskip("trimesh", reason="trimesh not installed, test skipped")
    cpt_mesh = Mesh.from_list_of_faces([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]])
    trimesh_mesh = cpt_mesh.export_to_trimesh()
    assert trimesh_mesh.faces.shape == (2, 3)


def test_compressed_stl():
    pytest.importorskip("trimesh", reason="trimesh not installed, test skipped")
    path = os.path.join(
        os.path.dirname(__file__),
        "mesh_files_examples/viking_ship.stl.xz"
    )
    with lzma.open(path, 'r') as f:
        mesh = load_mesh(f, file_format="stl", backend="trimesh")
    assert mesh.nb_faces == 2346
