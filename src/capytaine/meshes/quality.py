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
from itertools import chain

import numpy as np

LOG = logging.getLogger(__name__)
SHOWED_PYVISTA_WARNING = False


def check_mesh_quality(mesh, *, tol=1e-8):
    """
    Perform a set of geometric and metric quality checks on mesh data.

    Checks performed:
    - Non-coplanar faces
    - Non-convex faces
    - Aspect ratio via PyVista (if available)
    """
    non_coplanar = indices_of_non_coplanar_faces(mesh.vertices, mesh._faces)
    if non_coplanar:
        LOG.warning(f"{len(non_coplanar)} non-coplanar faces detected in {mesh}.")

    non_convex = indices_of_non_convex_faces(mesh.vertices, mesh._faces)
    if non_convex:
        LOG.warning(f"{len(non_convex)} non-convex faces detected in {mesh}.")

    try:
        pv_mesh = mesh.export_to_pyvista()

        try:
            arq = pv_mesh.cell_quality("aspect_ratio").cell_data.get("aspect_ratio")
        except AttributeError:  # Older version of PyVista
            arq = pv_mesh.compute_cell_quality("aspect_ratio").cell_data.get("CellQuality")
        if arq is not None:
            ratio_ok = np.sum(arq < 5) / len(arq)
            if ratio_ok < 0.9:
                LOG.info(f"Aspect Ratio of {mesh}:")
                LOG.info(
                    f"  Min: {np.min(arq):.3f} | Max: {np.max(arq):.3f} | Mean: {np.mean(arq):.3f}"
                )
                LOG.info(f"  Elements with AR < 5: {ratio_ok*100:.1f}%")
                LOG.warning(
                    "Low quality: more than 10% of elements have aspect ratio higher than 5."
                )
    except ImportError:
        global SHOWED_PYVISTA_WARNING
        if LOG.isEnabledFor(logging.INFO) and not SHOWED_PYVISTA_WARNING:
            LOG.info("PyVista not installed, skipping aspect ratio check.")
            SHOWED_PYVISTA_WARNING = True


def extract_face_vertices(vertices, face):
    return vertices[face]


def is_non_coplanar(vertices):
    a, b, c, d = vertices
    normal = np.cross(b - a, c - a)
    deviation = np.abs(np.dot(normal, d - a))
    return deviation > 1e-8


def indices_of_non_coplanar_faces(vertices, faces):
    """
    Identify the indices of quadrilateral faces that are not coplanar.

    Parameters
    ----------
    vertices : np.ndarray
        Array of vertex coordinates (n_vertices, 3).
    faces : np.ndarray
        Array of face indices (n_faces, 4) or (n_faces, 3).

    Returns
    -------
    list[int]
        List of indices of non-coplanar quadrilateral faces.
    """
    indices = []
    for i, face in enumerate(faces):
        if len(face) != 4:
            continue  # skip triangles or degenerate quads
        verts = extract_face_vertices(vertices, face)
        if is_non_coplanar(verts):
            indices.append(i)
    return indices


def is_face_convex(vertices):
    a, b, c, d = vertices
    edges = [b - a, c - b, d - c, a - d]
    normals = [np.cross(edges[i], edges[(i + 1) % 4]) for i in range(4)]
    dot_signs = [np.dot(normals[0], n) for n in normals[1:]]
    return all(s >= -1e-10 for s in dot_signs)


def indices_of_non_convex_faces(vertices, faces):
    """
    Identify indices of quadrilateral faces in the mesh that are not convex.

    Parameters
    ----------
    mesh : Mesh
        The input mesh containing faces and vertices.

    Returns
    -------
    list[int]
        List of indices of non-convex quadrilateral faces.
    """
    indices = []
    for i, face in enumerate(faces):
        if len(face) != 4:
            continue

        verts = extract_face_vertices(vertices, face)

        if not is_face_convex(verts):
            indices.append(i)

    return indices


def _is_valid(vertices, faces):
    """
    Check that every face index is valid for the given vertices.
    Now: an empty face-list is considered valid (as long as vertices exist).
    """
    # If you have no vertices at all, only accept if also no faces
    if len(vertices) == 0:
        return len(faces) == 0

    # If you have vertices but zero faces → that's fine
    if len(faces) == 0:
        return True

    # Otherwise, flatten all face‐indices and check bounds
    all_idx = list(chain.from_iterable(faces))
    if not all_idx:
        return True

    if min(all_idx) < 0 or max(all_idx) >= len(vertices):
        return False

    return True
