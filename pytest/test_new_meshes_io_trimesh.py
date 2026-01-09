
import pytest

from capytaine.new_meshes.io import load_mesh
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh

def test_load_directly_from_trimesh(tmp_path):
    """
    Test loading a simple triangle mesh from an .obj file using trimesh.
    """
    pytest.importorskip("trimesh")

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
    pytest.importorskip("trimesh")

    import trimesh
    trimesh_mesh = trimesh.Trimesh(
        vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], faces=[[0, 1, 2]]
    )

    test_file = tmp_path / "dummy.glb"
    trimesh_mesh.export(test_file)

    mesh = load_mesh(test_file, "glb")
    assert mesh.nb_faces == 1
    assert mesh.nb_vertices == 3
