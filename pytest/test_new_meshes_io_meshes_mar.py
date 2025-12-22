"""Tests for the mesh import of mesh from DIODORE-MARE file format"""

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
