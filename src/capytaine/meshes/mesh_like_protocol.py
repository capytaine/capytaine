from typing import Tuple, Protocol, runtime_checkable
from numpy.typing import ArrayLike


@runtime_checkable
class MeshLike(Protocol):
    """Minimal API that a class describing a mesh should implement to be
    usable with the rest of Capytaine.

    The goal is two-fold:
       1. Use at runtime to identify a mesh for functions that behaves
       differently depending on the type of the input (e.g. Delhommeau().evaluate).
       2. Use as documentation for third-party mesh implementation.

    In the future, it could also be used for static typing.
    """
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

    def join_meshes(*meshes, return_mask):
        ...

    def with_normal_vector_going_down(self, **kwargs):
        ...
