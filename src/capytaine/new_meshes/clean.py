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

import logging
from typing import List, Tuple, Dict

import numpy as np
from scipy.spatial import cKDTree

LOG = logging.getLogger(__name__)


def clean_mesh(
    vertices: np.ndarray,
    faces: List[List[int]],
    faces_metadata: Dict[str, np.ndarray],
    max_iter: int = 5,
    tol: float = 1e-8
) -> Tuple[np.ndarray, List[List[int]]]:
    """Iteratively clean a mesh by applying geometric simplifications.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of the input mesh.
    faces : list of list of int
        Face connectivity describing the mesh panels.
    faces_metadata: Dict[str, np.ndarray]
        Some arrays with the same first dimension (should be the number of faces)
        storing some fields defined on all the faces of the mesh.
    max_iter : int, default=5
        Maximum number of cleaning iterations to perform.
    tol : float, default=1e-8
        Tolerance used when merging near-duplicate vertices.

    Returns
    -------
    tuple[numpy.ndarray, list of list of int]
        The cleaned vertex array and associated face connectivity.
    """
    for _ in range(max_iter):
        nb_vertices_before = len(vertices)
        nb_faces_before = len(faces)

        vertices, faces, faces_metadata = clean_mesh_once(vertices, faces, faces_metadata, tol=tol)

        if len(vertices) == nb_vertices_before and len(faces) == nb_faces_before:
            break

    return vertices, faces, faces_metadata


def clean_mesh_once(
    vertices: np.ndarray,
    faces: List[List[int]],
    faces_metadata: Dict[str, np.ndarray],
    tol: float = 1e-10
) -> Tuple[np.ndarray, List[List[int]]]:
    """Run a single cleaning pass on the mesh data.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates describing the mesh geometry.
    faces : list of list of int
        Face connectivity with indices referencing ``vertices``.
    faces_metadata: Dict[str, np.ndarray]
        Some arrays with the same first dimension (should be the number of faces)
        storing some fields defined on all the faces of the mesh.
    tol : float, default=1e-10
        Tolerance for considering vertices as duplicates.

    Returns
    -------
    tuple[numpy.ndarray, list of list of int]
        Updated vertices and faces after the cleaning step.

    Raises
    ------
    ValueError
        If an unsupported face configuration is encountered.
    """
    # 1) merge almost‐duplicate vertices
    vertices, faces = merge_near_duplicate_vertices(vertices, faces, tol=tol)

    # 2) remove duplicate vertices indices in faces
    # and check that all faces have 3 or 4 unique vertices
    new_faces = []
    degenerate_faces_indices = []

    for i_face, face in enumerate(faces):
        seen = set()
        uniq = []
        for vi in face:
            if vi not in seen:
                seen.add(vi)
                uniq.append(vi)

        if len(uniq) in (3, 4):
            new_faces.append(uniq)
        elif len(uniq) < 3:
            degenerate_faces_indices.append(i_face)
        else:
            raise ValueError(
                f"Face with {len(uniq)} unique vertices: only 3 or 4 supported."
            )

    if len(degenerate_faces_indices) > 0:
        LOG.warning(
            f"Dropping {len(degenerate_faces_indices)} degenerate faces with <3 vertices: "
            f"{[faces[i] for i in degenerate_faces_indices[:5]]}{' ...' if len(degenerate_faces_indices) > 5 else ''}"
        )
        faces_metadata = {k: np.delete(faces_metadata[k], degenerate_faces_indices, axis=0) for k in faces_metadata}

    # 3) continue cleaning pipeline, all functions must accept List-of-lists too
    vertices, faces = remove_duplicate_vertices(vertices, new_faces)
    faces, faces_metadata = remove_duplicate_faces(faces, faces_metadata)
    vertices, faces = remove_unused_vertices(vertices, faces)
    faces, faces_metadata = remove_small_faces(vertices, faces, faces_metadata, tol=tol)
    vertices, faces = remove_unused_vertices(vertices, faces)

    return vertices, faces, faces_metadata


def merge_near_duplicate_vertices(
    vertices: np.ndarray, faces: List[List[int]], tol: float = 1e-8
) -> Tuple[np.ndarray, List[List[int]]]:
    """Merge vertices that are closer than a tolerance.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of shape ``(n, 3)``.
    faces : list of list of int
        Face connectivity referencing the ``vertices`` array.
    tol : float, default=1e-8
        Distance threshold below which vertices are considered duplicates.

    Returns
    -------
    tuple[numpy.ndarray, list of list of int]
        Deduplicated vertices and remapped faces.
    """
    if len(vertices) == 0:
        return vertices, faces

    tree = cKDTree(vertices)
    groups = tree.query_ball_tree(tree, r=tol)

    representative = {}
    new_vertices = []
    for i, group in enumerate(groups):
        rep = min(group)
        if rep not in representative:
            representative[rep] = len(new_vertices)
            new_vertices.append(vertices[rep])
        representative[i] = representative[rep]

    faces = [[representative[idx] for idx in face] for face in faces]
    new_vertices = np.array(new_vertices)
    return new_vertices, faces


def remove_duplicate_vertices(
    vertices: np.ndarray, faces: List[List[int]]
) -> Tuple[np.ndarray, List[List[int]]]:
    """Remove exactly repeated vertices and remap faces accordingly.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of shape ``(n, 3)``.
    faces : list of list of int
        Face connectivity using indices into ``vertices``.

    Returns
    -------
    tuple[numpy.ndarray, list of list of int]
        Unique vertices and faces with updated indices.
    """
    unique_vertices = []
    vertices_map = {}
    for vertex in vertices:
        vertex_tuple = tuple(vertex)
        if vertex_tuple not in vertices_map:
            vertices_map[vertex_tuple] = len(unique_vertices)
            unique_vertices.append(vertex)
    new_faces = [[vertices_map[tuple(vertices[i])] for i in face] for face in faces]
    new_vertices = np.array(unique_vertices)

    return new_vertices, new_faces


def remove_duplicate_faces(faces, faces_metadata):
    """Eliminate duplicate faces while preserving order.

    Parameters
    ----------
    faces : list of list of int
        Face connectivity to deduplicate.
    faces_metadata: Dict[str, np.ndarray]
        Fields associated to faces

    Returns
    -------
    list of list of int
        Face connectivity with duplicates removed.
    Dict[str, np.ndarray]
        Updated metadata
    """
    unique_faces = []
    face_set = set()
    deduplicated_faces_indices = []
    for i_face, face in enumerate(faces):
        face_tuple = tuple(sorted(face))
        if face_tuple not in face_set:
            face_set.add(face_tuple)
            unique_faces.append(face)
        else:
            deduplicated_faces_indices.append(i_face)

    faces_metadata = {k: np.delete(faces_metadata[k], deduplicated_faces_indices, axis=0) for k in faces_metadata}

    return unique_faces, faces_metadata


def remove_unused_vertices(
    vertices: np.ndarray, faces: List[List[int]]
) -> Tuple[np.ndarray, List[List[int]]]:
    """Remove vertices that are not referenced by any face.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of shape ``(n, 3)``.
    faces : list of list of int
        Face connectivity using indices into ``vertices``.

    Returns
    -------
    tuple[numpy.ndarray, list of list of int]
        Reduced vertex array and corresponding face connectivity.
    """
    used = sorted({i for face in faces for i in face})
    remap = {old: new for new, old in enumerate(used)}
    new_vs = vertices[used]
    new_fs = [[remap[i] for i in face] for face in faces]
    return new_vs, new_fs


def remove_small_faces(
    vertices: np.ndarray,
    faces: List[List[int]],
    faces_metadata: Dict[str, np.ndarray],
    tol: float = 1e-8
):
    """Remove faces whose area falls below a tolerance.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates used to evaluate surface area.
    faces : list of list of int
        Face connectivity referencing ``vertices``.
    faces_metadata: Dict[str, np.ndarray]
        Fields associated to faces
    tol : float, default=1e-8
        Minimum allowable face area.

    Returns
    -------
    list of list of int
        Faces that exceed the area threshold.
    """

    def face_area(face):
        v = vertices[face]
        if len(face) == 4:
            a1 = 0.5 * np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0]))
            a2 = 0.5 * np.linalg.norm(np.cross(v[2] - v[0], v[3] - v[0]))
            return a1 + a2
        elif len(face) == 3:
            return 0.5 * np.linalg.norm(np.cross(v[1] - v[0], v[2] - v[0]))
        return 0.0

    areas = np.array([face_area(face) for face in faces])
    mask = areas > tol
    faces = [face for face, keep in zip(faces, mask) if keep]
    faces_metadata = {k: faces_metadata[k][mask, ...] for k in faces_metadata}

    return faces, faces_metadata
