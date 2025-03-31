from typing import Tuple, Protocol, runtime_checkable
from numpy.typing import ArrayLike


@runtime_checkable
class MeshLike(Protocol):
    """Minimal API that a class describing a mesh should implement to be
    usable with the rest of Capytaine."""
    vertices: ArrayLike
    faces: ArrayLike
    nb_vertices: int
    nb_faces: int
    faces_centers: ArrayLike
    faces_normals: ArrayLike
    faces_areas: ArrayLike
    faces_radiuses: ArrayLike
    quadrature_points: Tuple[ArrayLike, ArrayLike]

    def __short_str__(self) -> str:
        ...

    def extract_faces(self, faces_id):
        ...
