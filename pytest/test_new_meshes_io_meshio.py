import os
import pytest

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

    mesh = load_mesh(test_file, "obj")
    assert mesh.nb_faces == 3
    assert mesh.nb_vertices == 6


def test_MED_file():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
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
