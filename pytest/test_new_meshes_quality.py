import capytaine.new_meshes as meshes
from capytaine.new_meshes.quality import (
    indices_of_non_convex_faces,
    indices_of_non_coplanar_faces,
)


def test_coplanar_quads():
    """All quadrilaterals are coplanar → no indices returned."""
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ]
    )
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == []


def test_one_non_coplanar_quad():
    """First quadrilateral is non-coplanar → index [0]."""
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ]
    )
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == [0]


def test_two_non_coplanar_quads():
    """Both quadrilaterals are non-coplanar → indices [0, 1]."""
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.1], [1.0, 1.0, 0.0]],
        ]
    )
    assert indices_of_non_coplanar_faces(mesh.vertices, mesh.faces) == [0, 1]


def test_convex_quads_are_valid():
    """Convex quadrilaterals should not be flagged as non-convex."""
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ]
    )
    assert indices_of_non_convex_faces(mesh.vertices, mesh.faces) == []


def test_non_convex_quad_is_detected():
    """A quadrilateral with crossed diagonals (vertex order issue) should be flagged."""
    mesh = meshes.Mesh.from_list_of_faces(
        [
            [
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],  # swapped vertices → concave / self-intersecting
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ]
    )
    assert indices_of_non_convex_faces(mesh.vertices, mesh.faces) == [0]
