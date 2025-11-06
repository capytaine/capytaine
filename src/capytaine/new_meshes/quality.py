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
from itertools import chain

import numpy as np

LOG = logging.getLogger(__name__)


def check_mesh_quality(vertices, faces, *, tol=1e-8):
    """
    Perform a set of geometric and metric quality checks on mesh data.

    Checks performed:
    - Duplicate vertices
    - Non-coplanar faces
    - Non-convex faces
    - Aspect ratio via PyVista (if available)

    Raises
    ------
    ValueError
        If mesh invalid data are provided.
    """
    non_coplanar = indices_of_non_coplanar_faces(vertices, faces)
    if non_coplanar:
        LOG.warning(f"{len(non_coplanar)} non-coplanar faces detected.")

    non_convex = indices_of_non_convex_faces(vertices, faces)
    if non_convex:
        LOG.warning(f"{len(non_convex)} non-convex faces detected.")

    warn_superimposed_faces(vertices, faces, tol=tol)

    try:
        pv_mesh = mesh_to_pyvista(vertices, faces)

        arq = pv_mesh.compute_cell_quality("aspect_ratio").cell_data.get("CellQuality")
        if arq is not None:
            ratio_ok = np.sum(arq < 5) / len(arq)
            if ratio_ok < 0.9:
                LOG.info("Aspect Ratio:")
                LOG.info(
                    f"  Min: {np.min(arq):.3f} | Max: {np.max(arq):.3f} | Mean: {np.mean(arq):.3f}"
                )
                LOG.info(f"  Elements with AR < 5: {ratio_ok*100:.1f}%")
                LOG.warning(
                    "Low quality: more than 10%% of elements have aspect ratio higher than 5."
                )
    except ModuleNotFoundError:
        LOG.info("PyVista not installed, skipping aspect ratio check.")


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


def mesh_to_pyvista(
    vertices: np.ndarray, faces: List[List[int]]
) -> "pv.UnstructuredGrid":
    """
    Build a PyVista UnstructuredGrid from a list of faces (triangles or quads).
    """
    import pyvista as pv
    # flatten into the VTK cell‐array format: [n0, i0, i1, ..., in-1, n1, j0, j1, ...]
    flat_cells = []
    cell_types = []
    for face in faces:
        n = len(face)
        flat_cells.append(n)
        flat_cells.extend(face)
        if n == 3:
            cell_types.append(pv.CellType.TRIANGLE)
        elif n == 4:
            cell_types.append(pv.CellType.QUAD)
        else:
            # if you ever have ngons, you can map them as POLYGON:
            cell_types.append(pv.CellType.POLYGON)

    cells_array = np.array(flat_cells, dtype=np.int64)
    cell_types = np.array(cell_types, dtype=np.uint8)

    return pv.UnstructuredGrid(cells_array, cell_types, vertices.astype(np.float32))


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


def warn_superimposed_faces(
    vertices: np.ndarray, faces: List[List[int]], *, tol: float = 1e-8
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

    keyed_faces: Dict[Tuple, Tuple[int, np.ndarray]] = {}
    duplicates: List[Tuple[int, int, float]] = []
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
