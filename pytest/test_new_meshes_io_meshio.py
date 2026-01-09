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

from capytaine.new_meshes.meshes import Mesh
from capytaine.new_meshes.io import load_mesh

def test_load_directly_from_meshio():
    pytest.importorskip("meshio")

    import meshio
    points = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
    ]
    cells = [("triangle", [[0, 1, 2], [1, 3, 2]]), ("quad", [[1, 4, 5, 3]])]
    meshio_mesh = meshio.Mesh(points, cells)

    mesh = load_mesh(meshio_mesh)
    assert mesh.nb_faces == 3
    assert mesh.nb_vertices == 6


def test_load_from_meshio(tmp_path):
    pytest.importorskip("meshio")

    import meshio
    points = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
    ]
    cells = [("triangle", [[0, 1, 2], [1, 3, 2]]), ("quad", [[1, 4, 5, 3]])]
    meshio_mesh = meshio.Mesh(points, cells)

    mesh = load_mesh(meshio_mesh)

    test_file = tmp_path / "dummy.obj"
    meshio.write(test_file, meshio_mesh)

    mesh = load_mesh(test_file, "obj", backend="meshio")
    assert mesh.nb_faces == 3
    assert mesh.nb_vertices == 6


def test_MED_file():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    pytest.importorskip("h5py", reason="h5py not installed, test skipped")
    mesh = load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/barge.med"))
    assert mesh.nb_faces == 187

def test_MSH2_path():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    mesh = load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder2.msh"), file_format="gmsh")
    assert mesh.nb_faces == 64

def test_MSH2_file():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    with open(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder2.msh"), 'r') as f:
        mesh = load_mesh(f, file_format="gmsh")
    assert mesh.nb_faces == 64

def test_MSH4_path():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    mesh = load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder4.msh"), file_format="gmsh")
    assert mesh.nb_faces == 64

def test_MSH4_file():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    with open(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder4.msh"), 'r') as f:
        mesh = load_mesh(f, file_format="gmsh")
    assert mesh.nb_faces == 64

def test_export_to_meshio():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    cpt_mesh = Mesh.from_list_of_faces([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]])
    meshio_mesh = cpt_mesh.export_to_meshio()
    assert len(meshio_mesh.cells[0]) == 1

def test_compressed_stl():
    path = os.path.join(
        os.path.dirname(__file__),
        "mesh_files_examples/viking_ship.stl.xz"
    )
    with lzma.open(path, 'rb') as f:
        mesh = load_mesh(f, file_format="stl", backend="meshio")
    assert mesh.nb_faces == 2346
