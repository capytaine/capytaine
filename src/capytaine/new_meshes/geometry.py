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

from typing import List
from functools import reduce
from itertools import chain

import numpy as np
from numpy.typing import NDArray


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


def compute_gauss_legendre_2_quadrature(vertices, faces):
    # Parameters of Gauss-Legendre 2 quadrature scheme
    local_points = np.array([(+1/np.sqrt(3), +1/np.sqrt(3)),
             (+1/np.sqrt(3), -1/np.sqrt(3)),
             (-1/np.sqrt(3), +1/np.sqrt(3)),
             (-1/np.sqrt(3), -1/np.sqrt(3))])
    local_weights = np.array([1/4, 1/4, 1/4, 1/4])

    # Application to mesh
    faces = vertices[faces[:, :], :]
    nb_faces = faces.shape[0]
    nb_quad_points = len(local_weights)
    points = np.empty((nb_faces, nb_quad_points, 3))
    weights = np.empty((nb_faces, nb_quad_points))
    for i_face in range(nb_faces):
        for k_quad in range(nb_quad_points):
            xk, yk = local_points[k_quad, :]
            points[i_face, k_quad, :] = (
                      (1+xk)*(1+yk) * faces[i_face, 0, :]
                    + (1+xk)*(1-yk) * faces[i_face, 1, :]
                    + (1-xk)*(1-yk) * faces[i_face, 2, :]
                    + (1-xk)*(1+yk) * faces[i_face, 3, :]
                    )/4
            dxidx = ((1+yk)*faces[i_face, 0, :] + (1-yk)*faces[i_face, 1, :]
                     - (1-yk)*faces[i_face, 2, :] - (1+yk)*faces[i_face, 3, :])/4
            dxidy = ((1+xk)*faces[i_face, 0, :] - (1+xk)*faces[i_face, 1, :]
                     - (1-xk)*faces[i_face, 2, :] + (1-xk)*faces[i_face, 3, :])/4
            detJ = np.linalg.norm(np.cross(dxidx, dxidy))
            weights[i_face, k_quad] = local_weights[k_quad] * 4 * detJ

    return points, weights


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

def faces_in_group(faces: NDArray[np.integer], group: NDArray[np.integer]) -> NDArray[np.bool_]:
    """Identification of faces with vertices within group.

    Parameters
    ----------
    faces : NDArray[np.integer]
        Mesh faces. Expecting a numpy array of shape N_faces x N_vertices_per_face.
    group : NDArray[np.integer]
        Group of connected vertices

    Returns
    -------
    NDArray[np.bool]
        Mask of faces containing vertices from the group
    """
    return np.any(np.isin(faces, group), axis=1)

def clustering(faces: NDArray[np.integer]) -> List[NDArray[np.integer]]:
    """Clustering of vertices per connected faces.

    Parameters
    ----------
    faces : NDArray[np.integer]
        Mesh faces. Expecting a numpy array of shape N_faces x N_vertices_per_face.

    Returns
    -------
    list[NDArray[np.integer]]
        Groups of connected vertices.
    """
    vert_groups: list[NDArray[np.integer]] = []
    mask = np.ones(faces.shape[0], dtype=bool)
    while np.any(mask):
        # Consider faces whose vertices are not already identified in a group.
        # Start new group by considering first face
        remaining_faces = faces[mask]
        group = remaining_faces[0]
        rem_mask = np.ones(remaining_faces.shape[0], dtype=bool)
        # Iterative update of vertices group. Output final result to frozenset
        while not np.allclose(new:=faces_in_group(remaining_faces, group), rem_mask):
            group = np.unique(remaining_faces[new])
            rem_mask = new
        else:
            group = np.unique(remaining_faces[new])
        vert_groups.append(group)
        # Identify faces that have no vertices in current groups
        mask = ~reduce(np.logical_or, [faces_in_group(faces, group) for group in vert_groups])
    return vert_groups


def connected_components(mesh):
    """Returns a list of meshes that each corresponds to the a connected component in the original mesh.
    Assumes the mesh is mostly conformal without duplicate vertices.
    """
    # Get connected vertices
    vertices_components = clustering(mesh.faces)
    # Verification
    if sum(len(group) for group in vertices_components) != len(set(chain.from_iterable(vertices_components))):
        raise ValueError("Error in connected components clustering. Some elements are duplicated")
    # The components are found. The rest is just about retrieving the faces in each components.
    faces_components = [np.argwhere(faces_in_group(mesh.faces, group)) for group in vertices_components]
    components = [mesh.extract_faces(f) for f in faces_components]
    return components


def connected_components_of_waterline(mesh, z=0.0):
    if np.any(mesh.vertices[:, 2] > z + 1e-8):
        mesh = mesh.immersed_part(free_surface=z)
    fs_vertices_indices = np.where(np.isclose(mesh.vertices[:, 2], z))[0]
    fs_faces_indices = np.where(np.any(np.isin(mesh.faces, fs_vertices_indices), axis=1))[0]
    crown_mesh = mesh.extract_faces(fs_faces_indices)
    return connected_components(crown_mesh)
