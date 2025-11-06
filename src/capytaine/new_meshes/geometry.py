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

import numpy as np


def get_vertices_face(face, vertices):
    if len(face) == 4 and face[2] != face[3]:
        return (
            vertices[face[0]],
            vertices[face[1]],
            vertices[face[2]],
            vertices[face[3]],
        )
    else:
        return vertices[face[0]], vertices[face[1]], vertices[face[2]]


def compute_faces_normals(vertices, faces):
    normals = []
    for face in faces:
        if len(face) == 4 and face[2] != face[3]:
            normal = _quad_normal(vertices, face[0], face[1], face[2], face[3])
        else:
            normal = _triangle_normal(vertices, face[0], face[1], face[2])
        normals.append(normal)
    return np.array(normals)


def compute_faces_areas(vertices, faces):
    areas = []
    for face in faces:
        verts = get_vertices_face(face, vertices)
        if len(verts) == 4:
            a, b, c, d = verts
            area1 = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
            area2 = 0.5 * np.linalg.norm(np.cross(c - a, d - a))
            areas.append(area1 + area2)
        else:
            a, b, c = verts
            areas.append(0.5 * np.linalg.norm(np.cross(b - a, c - a)))
    return np.array(areas)


def compute_faces_centers(vertices, faces):
    centers = []
    for face in faces:
        verts = get_vertices_face(face, vertices)
        if len(verts) == 4:
            a, b, c, d = verts
            area1 = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
            area2 = 0.5 * np.linalg.norm(np.cross(c - a, d - a))
            c1 = (a + b + c) / 3
            c2 = (a + c + d) / 3
            center = (c1 * area1 + c2 * area2) / (area1 + area2)
        else:
            a, b, c = verts
            center = (a + b + c) / 3
        centers.append(center)
    return np.array(centers)


def compute_faces_radii(vertices, faces):
    centers = compute_faces_centers(vertices, faces)
    distances = []
    for face, center in zip(faces, centers):
        d = compute_distance_between_points(vertices[face[0]], center)
        distances.append(d)
    return np.array(distances)


def _triangle_normal(vertices, v0_idx, v1_idx, v2_idx):
    """
    Compute normal vector of a triangle face.

    Parameters
    ----------
    vertices : ndarray
        Vertex coordinate array.
    v0_idx, v1_idx, v2_idx : int
        Indices of triangle vertices.

    Returns
    -------
    np.ndarray
        Normalized normal vector (3,)
    """
    v0, v1, v2 = vertices[v0_idx], vertices[v1_idx], vertices[v2_idx]
    normal = np.cross(v1 - v0, v2 - v0)
    return normal / np.linalg.norm(normal)


def _quad_normal(vertices, v0_idx, v1_idx, v2_idx, v3_idx):
    """
    Compute normal vector of a quadrilateral face via diagonals.

    Parameters
    ----------
    vertices : ndarray
        Vertex coordinate array.
    v0_idx, v1_idx, v2_idx, v3_idx : int
        Indices of quad vertices.

    Returns
    -------
    np.ndarray
        Normalized normal vector (3,)
    """
    v0, v1, v2, v3 = (
        vertices[v0_idx],
        vertices[v1_idx],
        vertices[v2_idx],
        vertices[v3_idx],
    )
    ac = v2 - v0
    bd = v3 - v1
    normal = np.cross(ac, bd)
    return normal / np.linalg.norm(normal)


def compute_distance_between_points(a, b):
    """
    Compute Euclidean distance between two points in n-dimensional space.

    Parameters
    ----------
    a, b : array_like
        Coordinate arrays (length 3 or more).

    Returns
    -------
    float
        Euclidean distance.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    return np.linalg.norm(b - a)
