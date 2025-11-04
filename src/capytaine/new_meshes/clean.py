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
from typing import List

import numpy as np
from scipy.spatial import cKDTree

LOG = logging.getLogger(__name__)


def clean_mesh(
    vertices: np.ndarray, faces: List[List[int]], max_iter: int = 5, tol: float = 1e-8
) -> tuple[np.ndarray, List[List[int]]]:
    """Iteratively clean a mesh by applying geometric simplifications.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of the input mesh.
    faces : list of list of int
        Face connectivity describing the mesh panels.
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

        vertices, faces = clean_mesh_once(vertices, faces, tol=tol)

        if len(vertices) == nb_vertices_before and len(faces) == nb_faces_before:
            break

    return vertices, faces


def clean_mesh_once(
    vertices: np.ndarray, faces: List[List[int]], tol: float = 1e-10
) -> tuple[np.ndarray, List[List[int]]]:
    """Run a single cleaning pass on the mesh data.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates describing the mesh geometry.
    faces : list of list of int
        Face connectivity with indices referencing ``vertices``.
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

    # 2) collapse degenerate quads → tris, drop <3‐pt faces
    new_faces = []
    degenerate_faces = []

    for face in faces:
        seen = set()
        uniq = []
        for vi in face:
            if vi not in seen:
                seen.add(vi)
                uniq.append(vi)

        if len(uniq) in (3, 4):
            new_faces.append(uniq)
        elif len(uniq) < 3:
            degenerate_faces.append(uniq)
        else:
            raise ValueError(
                f"Face with {len(uniq)} unique vertices: only 3 or 4 supported."
            )

    if degenerate_faces:
        LOG.warning(
            f"Dropping {len(degenerate_faces)} degenerate faces with <3 vertices: "
            f"{degenerate_faces[:5]}{' ...' if len(degenerate_faces) > 5 else ''}"
        )

    warn_superimposed_faces(vertices, new_faces, tol=tol)

    # 3) continue cleaning pipeline, all functions must accept List-of-lists too
    vertices, faces = remove_duplicate_vertices(vertices, new_faces)
    faces = remove_duplicate_faces(faces)
    vertices, faces = remove_unused_vertices(vertices, faces)
    faces = remove_small_faces(vertices, faces, tol=tol)
    vertices, faces = remove_unused_vertices(vertices, faces)

    return vertices, faces


def merge_near_duplicate_vertices(
    vertices: np.ndarray, faces: List[List[int]], tol: float = 1e-8
) -> tuple[np.ndarray, List[List[int]]]:
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
) -> tuple[np.ndarray, List[List[int]]]:
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


def remove_duplicate_faces(faces: List[List[int]]) -> List[List[int]]:
    """Eliminate duplicate faces while preserving order.

    Parameters
    ----------
    faces : list of list of int
        Face connectivity to deduplicate.

    Returns
    -------
    list of list of int
        Face connectivity with duplicates removed.
    """
    unique_faces = []
    face_set = set()
    for face in faces:
        face_tuple = tuple(sorted(face))
        if face_tuple not in face_set:
            face_set.add(face_tuple)
            unique_faces.append(face)
    return unique_faces


def warn_superimposed_faces(
    vertices: np.ndarray, faces: List[List[int]], tol: float = 1e-8
):
    """Emit a warning when panels are duplicated within a tolerance.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates of the mesh.
    faces : list of list of int
        Face connectivity referencing ``vertices``.
    tol : float, default=1e-8
        Tolerance used to determine whether two faces coincide.
    """

    if not faces or len(faces) < 2:
        return

    keyed_faces: dict[tuple, tuple[int, np.ndarray]] = {}
    duplicates: list[tuple[int, int, float]] = []
    scale = 1.0 / max(tol, 1e-12)

    for idx, face in enumerate(faces):
        coords = vertices[face]
        quantized_points = [
            tuple(np.round(pt * scale).astype(np.int64)) for pt in coords
        ]
        quantized = tuple(sorted(quantized_points))
        sorted_coords = np.array(sorted(coords.tolist()))
        key = (len(face), quantized)
        if key in keyed_faces:
            base_idx, base_coords = keyed_faces[key]
            max_dev = float(
                np.max(np.linalg.norm(sorted_coords - base_coords, axis=1))
            )
            duplicates.append((base_idx, idx, max_dev))
        else:
            keyed_faces[key] = (idx, sorted_coords)

    if duplicates:
        sample = duplicates[:3]
        LOG.warning(
            "Detected %d panels that are superimposed or nearly identical (tol=%.1e). "
            "Example index pairs (i, j, max_dev): %s",
            len(duplicates),
            tol,
            sample,
        )


def remove_unused_vertices(
    vertices: np.ndarray, faces: List[List[int]]
) -> tuple[np.ndarray, List[List[int]]]:
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
    vertices: np.ndarray, faces: List[int], tol: float = 1e-8
) -> List[int]:
    """Remove faces whose area falls below a tolerance.

    Parameters
    ----------
    vertices : numpy.ndarray
        Vertex coordinates used to evaluate surface area.
    faces : list of int
        Face connectivity referencing ``vertices``.
    tol : float, default=1e-8
        Minimum allowable face area.

    Returns
    -------
    list of int
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

    return faces
