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

from __future__ import annotations
from typing import List, Set, Tuple

import numpy as np


def _get_intersection(
    v1: np.ndarray,
    v2: np.ndarray,
    normal: np.ndarray,
    origin: np.ndarray,
    tol: float = 1e-8,
) -> np.ndarray:
    """Intersect a line segment with a plane.

    Parameters
    ----------
    v1, v2 : numpy.ndarray
        Endpoints of the segment expressed as 3D vectors.
    normal : numpy.ndarray
        Plane normal; does not need to be unit length.
    origin : numpy.ndarray
        Point lying on the plane.
    tol : float, default=1e-8
        Tolerance used to detect parallel segments and clamp the intersection
        parameter to the segment bounds.

    Returns
    -------
    numpy.ndarray
        Intersection point located on the (possibly clamped) segment.

    Raises
    ------
    ValueError
        If the segment is parallel to the plane or the intersection lies
        outside the tolerated segment range.
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    n = np.asarray(normal, dtype=float)
    o = np.asarray(origin, dtype=float)

    u = v2 - v1
    denom = float(np.dot(n, u))
    if abs(denom) < tol:
        raise ValueError("Segment is parallel to the plane (no unique intersection).")

    t = float(np.dot(n, (o - v1)) / denom)

    if t < -tol or t > 1.0 + tol:
        raise ValueError(f"Intersection t={t:.6g} lies outside the segment [0,1].")

    t = max(0.0, min(1.0, t))
    return v1 + t * u


def compute_aspect_ratio(tri_pts: np.ndarray) -> float:
    """Compute the aspect ratio of a triangle.

    The aspect ratio is defined as the longest edge length divided by the
    altitude relative to that edge.

    Parameters
    ----------
    tri_pts : numpy.ndarray
        Triangle vertices arranged in a ``(3, 3)`` array.

    Returns
    -------
    float
        Aspect ratio ``L / h`` (values greater than or equal to 1).

    Raises
    ------
    ValueError
        If the triangle is degenerate and its area approaches zero.
    """
    tri_pts = np.asarray(tri_pts, dtype=float)
    if tri_pts.shape != (3, 3):
        raise ValueError("tri_pts must have shape (3,3).")

    edges = np.array(
        [
            np.linalg.norm(tri_pts[1] - tri_pts[0]),
            np.linalg.norm(tri_pts[2] - tri_pts[1]),
            np.linalg.norm(tri_pts[0] - tri_pts[2]),
        ],
        dtype=float,
    )
    L = float(np.max(edges))
    i = int(np.argmax(edges))

    A, B = tri_pts[i], tri_pts[(i + 1) % 3]
    C = tri_pts[(i + 2) % 3]
    area = 0.5 * np.linalg.norm(np.cross(B - A, C - A))
    if area <= 0.0:
        raise ValueError("Degenerate triangle: zero area.")
    h = 2.0 * area / L
    return L / h


def _signed_distances(
    verts: list[list[float]],
    face: list[int],
    normal: np.ndarray,
    origin: np.ndarray,
) -> np.ndarray:
    """Evaluate signed distances of face vertices to a clipping plane.

    Parameters
    ----------
    verts : list of list of float
        All mesh vertices.
    face : list of int
        Indices of the vertices forming the face.
    normal : numpy.ndarray
        Plane normal vector.
    origin : numpy.ndarray
        Point belonging to the plane.

    Returns
    -------
    numpy.ndarray
        Signed distances for the vertices belonging to ``face``.
    """
    pts = np.asarray([verts[i] for i in face], dtype=float)
    return (origin - pts) @ normal


def _compute_keep_sets(
    verts: list[list[float]],
    face: list[int],
    normal: np.ndarray,
    origin: np.ndarray,
    tol: float,
) -> tuple[list[int], set[int], np.ndarray]:
    """Split vertices of a face between kept and discarded sets.

    Parameters
    ----------
    verts : list of list of float
        All mesh vertices.
    face : list of int
        Indices forming the face currently being clipped.
    normal : numpy.ndarray
        Plane normal vector defining the clipping plane.
    origin : numpy.ndarray
        Point on the clipping plane.
    tol : float
        Tolerance used to consider vertices inside the kept half-space.

    Returns
    -------
    list of int
        Indices of vertices that remain after clipping.
    set of int
        Indices of vertices removed by the clipping plane.
    numpy.ndarray
        Signed distance of each face vertex to the plane.
    """
    s = _signed_distances(verts, face, normal, origin)
    mask = s >= -tol
    keep = [face[i] for i, m in enumerate(mask) if m]
    unkeep = set(face) - set(keep)
    return keep, unkeep, s


def _compute_edge_intersections(
    verts: list[list[float]],
    face: list[int],
    keep: list[int],
    normal: np.ndarray,
    origin: np.ndarray,
    tol: float,
) -> dict[tuple[int, int], int]:
    """Intersect face edges with the clipping plane.

    Parameters
    ----------
    verts : list of list of float
        Mutable list of mesh vertices; intersection points are appended here.
    face : list of int
        Indices of the vertices forming the face being processed.
    keep : list of int
        Vertices that remain inside the kept half-space.
    normal : numpy.ndarray
        Plane normal vector.
    origin : numpy.ndarray
        Point belonging to the plane.
    tol : float
        Tolerance for plane intersection checks.

    Returns
    -------
    dict[tuple[int, int], int]
        Mapping from directed edges to newly created vertex indices.
    """
    keep_set = set(keep)
    edges = list(zip(face, face[1:] + face[:1]))
    cut_edges = [(i, j) for (i, j) in edges if (i in keep_set) ^ (j in keep_set)]

    edge_inters: dict[tuple[int, int], int] = {}
    for i, j in cut_edges:
        ip = _get_intersection(
            np.array(verts[i]), np.array(verts[j]), normal, origin, tol
        )
        idx = len(verts)
        verts.append(ip.tolist())
        edge_inters[(i, j)] = idx
    return edge_inters


def _build_clipped_boundary(
    face: list[int],
    keep: list[int],
    edge_inters: dict[tuple[int, int], int],
) -> list[int]:
    """Assemble the boundary indices for a clipped polygon.

    Parameters
    ----------
    face : list of int
        Original face indices.
    keep : list of int
        Vertices that remain after clipping.
    edge_inters : dict[tuple[int, int], int]
        Mapping from edges to newly created intersection vertices.

    Returns
    -------
    list of int
        Ordered vertex indices describing the clipped boundary.
    """
    keep_set = set(keep)
    boundary: list[int] = []
    for i, j in zip(face, face[1:] + face[:1]):
        if i in keep_set:
            boundary.append(i)
        if (i, j) in edge_inters:
            boundary.append(edge_inters[(i, j)])
    return boundary


def clip_faces(
    vertices: np.ndarray,
    faces: list[list[int]],
    normal: np.ndarray,
    origin: np.ndarray,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, List[List[int]], np.ndarray]:
    """Clip faces of a mesh against a plane.

    The kept half-space is defined by ``(v - origin) · normal >= -tol``.

    Parameters
    ----------
    vertices : numpy.ndarray
        Input vertex positions of shape ``(n, 3)``.
    faces : list of list of int
        Face connectivity; triangles and quads are supported.
    normal : numpy.ndarray
        Normal vector of the clipping plane.
    origin : numpy.ndarray
        Point located on the plane.
    tol : float, default=1e-8
        Tolerance for classifying vertices relative to the plane.

    Returns
    -------
    np.ndarray of floats of shape (new_nb_vertices, 3)
        The pruned vertex array
    list of list of int
        The list of of length new_nb_faces with clipped faces
    np.ndarray of ints of shape (new_nb_faces,)
        For each new face, the index of the face it comes from in the input.
    """
    normal = np.asarray(normal, dtype=float)
    origin = np.asarray(origin, dtype=float)

    verts: List[List[float]] = vertices.astype(float).tolist()
    new_faces: List[Tuple[int, List[int]]] = []
    # A new face is a tuple storing the index of the parent face in the
    # original mesh and a list of vertices
    dropped_vs: Set[int] = set()

    for i_face, face in enumerate(faces):
        keep, unkeep, _ = _compute_keep_sets(verts, face, normal, origin, tol)
        dropped_vs.update(unkeep)

        if len(keep) == 0:
            continue  # fully outside
        if len(keep) == len(face):
            # Face fully inside → keep original (quad or triangle unchanged)
            new_faces.append((i_face, list(face)))
            continue

        edge_inters = _compute_edge_intersections(
            verts, face, keep, normal, origin, tol
        )
        boundary = _build_clipped_boundary(face, keep, edge_inters)

        if len(boundary) == 3:
            new_faces.append((i_face, boundary))
        elif len(boundary) == 4:
            # clipped quad → 2 triangles
            new_faces.append((i_face, [boundary[0], boundary[1], boundary[2]]))
            new_faces.append((i_face, [boundary[0], boundary[2], boundary[3]]))
        elif len(boundary) == 5:
            # pentagon → 1 triangle + 1 quad
            tri = [boundary[0], boundary[1], boundary[2]]
            quad = [boundary[0], boundary[2], boundary[3], boundary[4]]
            new_faces.append((i_face, tri))
            new_faces.append((i_face, quad))
        else:
            # fallback: fan triangulation
            for k in range(1, len(boundary) - 1):
                new_faces.append((i_face, [boundary[0], boundary[k], boundary[k + 1]]))

    if not new_faces:
        return np.empty((0, 3), dtype=float), [], np.empty((0,), dtype=int)

    used = {idx for (_, f) in new_faces for idx in f}
    dropped_vs -= used
    keep_vs = [i for i in range(len(verts)) if i not in dropped_vs]
    remap = {old: new for new, old in enumerate(keep_vs)}

    pruned_verts = np.asarray([verts[i] for i in keep_vs], dtype=float)
    pruned_faces = [[remap[i] for i in face] for (_, face) in new_faces]

    parent_of_face = np.array([i_parent_face for (i_parent_face, _) in new_faces])

    return pruned_verts, pruned_faces, parent_of_face
