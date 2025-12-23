import pytest

from capytaine.new_meshes import Mesh
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
