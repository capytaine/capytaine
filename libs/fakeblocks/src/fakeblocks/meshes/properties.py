"""Helper functions to compute some properties of the mesh.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

from functools import reduce
from itertools import chain
import numpy as np
from typing import List
from numpy.typing import NDArray


def compute_faces_properties(mesh):
    """Compute the faces properties of the mesh"""

    # faces_areas, faces_normals, faces_centers = mm.get_all_faces_properties(mesh._vertices, mesh._faces)
    nf = mesh.nb_faces

    # triangle_mask = _faces[:, 0] == _faces[:, -1]
    # nb_triangles = np.sum(triangle_mask)
    # quads_mask = np.invert(triangle_mask)
    # nb_quads = nf - nb_triangles

    faces_areas = np.zeros(nf, dtype=float)
    faces_normals = np.zeros((nf, 3), dtype=float)
    faces_centers = np.zeros((nf, 3), dtype=float)

    # Collectively dealing with triangles
    # triangles = _faces[triangle_mask]
    triangles_id = mesh.triangles_ids
    triangles = mesh._faces[triangles_id]

    triangles_normals = np.cross(mesh._vertices[triangles[:, 1]] - mesh._vertices[triangles[:, 0]],
                                 mesh._vertices[triangles[:, 2]] - mesh._vertices[triangles[:, 0]])
    triangles_normals_norm = np.linalg.norm(triangles_normals, axis=1)

    degenerate_triangle = np.abs(triangles_normals_norm) < 1e-12
    triangles_id = triangles_id[~degenerate_triangle]
    triangles_normals = triangles_normals[~degenerate_triangle, :]
    triangles_normals_norm = triangles_normals_norm[~degenerate_triangle]
    triangles = triangles[~degenerate_triangle, :]
    # Now, continue the computations without the degenerate triangles

    faces_normals[triangles_id] = triangles_normals / triangles_normals_norm[:, np.newaxis]
    faces_areas[triangles_id] = triangles_normals_norm / 2.
    faces_centers[triangles_id] = np.sum(mesh._vertices[triangles[:, :3]], axis=1) / 3.

    # Collectively dealing with quads
    quads_id = mesh.quadrangles_ids
    quads = mesh._faces[quads_id]
    # quads = _faces[quads_mask]

    quads_normals = np.cross(mesh._vertices[quads[:, 2]] - mesh._vertices[quads[:, 0]],
                             mesh._vertices[quads[:, 3]] - mesh._vertices[quads[:, 1]])

    quads_normals_norm = np.linalg.norm(quads_normals, axis=1)

    degenerate_quad = np.abs(quads_normals_norm) < 1e-12
    quads_id = quads_id[~degenerate_quad]
    quads_normals = quads_normals[~degenerate_quad]
    quads_normals_norm = quads_normals_norm[~degenerate_quad]
    quads = quads[~degenerate_quad, :]
    # Now, continue the computations without the degenerate quads

    faces_normals[quads_id] = quads_normals / quads_normals_norm[:, np.newaxis]

    a1 = np.linalg.norm(np.cross(mesh._vertices[quads[:, 1]] - mesh._vertices[quads[:, 0]],
                                 mesh._vertices[quads[:, 2]] - mesh._vertices[quads[:, 0]]), axis=1) * 0.5
    a2 = np.linalg.norm(np.cross(mesh._vertices[quads[:, 3]] - mesh._vertices[quads[:, 0]],
                                 mesh._vertices[quads[:, 2]] - mesh._vertices[quads[:, 0]]), axis=1) * 0.5
    faces_areas[quads_id] = a1 + a2

    c1 = np.sum(mesh._vertices[quads[:, :3]], axis=1) / 3.
    c2 = (np.sum(mesh._vertices[quads[:, 2:4]], axis=1) + mesh._vertices[quads[:, 0]]) / 3.

    faces_centers[quads_id] = (np.array(([a1, ] * 3)).T * c1 + np.array(([a2, ] * 3)).T * c2)
    faces_centers[quads_id] /= np.array(([faces_areas[quads_id], ] * 3)).T

    faces_radiuses = compute_radiuses(mesh, faces_centers)

    return {'faces_areas': faces_areas,
            'faces_normals': faces_normals,
            'faces_centers': faces_centers,
            'faces_radiuses': faces_radiuses,
            }


def compute_radiuses(mesh, faces_centers):
    """Compute the radiuses of the faces of the mesh.

    The radius is defined here as the maximal distance between the center
    of mass of a cell and one of its points."""

    # Coordinates of all the vertices grouped by face
    faces_vertices = mesh.vertices[mesh.faces, :]
    # faces_vertices.shape == (nb_faces, 4, 3)

    # Reorder the axes for array broadcasting below
    faces_vertices = np.moveaxis(faces_vertices, 0, 1)
    # faces_vertices.shape == (4, nb_faces, 3)

    # Get all the vectors between the center of faces and their vertices.
    radial_vector = faces_centers - faces_vertices
    # radial_vector.shape == (4, nb_faces, 3)

    # Keep the maximum length
    faces_radiuses = np.max(np.linalg.norm(radial_vector, axis=2), axis=0)
    # faces_radiuses.shape = (nb_faces)

    return faces_radiuses


def compute_connectivity(mesh):
    """Compute the connectivities of the mesh.

    It concerns further connectivity than simple faces/vertices connectivities. It computes the vertices / vertices, vertices / faces and faces / faces connectivities.

    Note
    ----
    * Note that if the mesh is not conformal, the algorithm may not perform correctly

    TODO: The computation of boundaries should be in the future computed with help of VTK
    """

    nv = mesh.nb_vertices
    nf = mesh.nb_faces

    mesh_closed = True

    # Building connectivities

    # Establishing v_v and v_f connectivities
    v_v = dict([(i, set()) for i in range(nv)])
    v_f = dict([(i, set()) for i in range(nv)])
    for (iface, face) in enumerate(mesh._faces):
        if face[0] == face[-1]:
            face_w = face[:3]
        else:
            face_w = face
        for (index, iV) in enumerate(face_w):
            v_f[iV].add(iface)
            v_v[face_w[index - 1]].add(iV)
            v_v[iV].add(face_w[index - 1])

    # Connectivity f_f
    boundary_edges = dict()

    f_f = dict([(i, set()) for i in range(nf)])
    for ivertex in range(nv):
        set1 = v_f[ivertex]
        for iadj_v in v_v[ivertex]:
            set2 = v_f[iadj_v]
            intersection = list(set1 & set2)
            if len(intersection) == 2:
                f_f[intersection[0]].add(intersection[1])
                f_f[intersection[1]].add(intersection[0])

            elif len(intersection) == 1:
                boundary_face = mesh._faces[intersection[0]]

                if boundary_face[0] == boundary_face[-1]:
                    boundary_face = boundary_face[:3]
                ids = np.where((boundary_face == ivertex) + (boundary_face == iadj_v))[0]

                if ids[1] != ids[0]+1:
                    i_v_orig, i_v_target = boundary_face[ids]
                else:
                    i_v_target, i_v_orig = boundary_face[ids]

                boundary_edges[i_v_orig] = i_v_target
            else:
                raise RuntimeError('Unexpected error while computing mesh connectivities')

    # Computing boundaries
    boundaries = list()
    # TODO: calculer des boundaries fermees et ouvertes (closed_boundaries et open_boundaries) et mettre dans dict
    while True:
        try:
            boundary = list()
            i_v0_init, i_v1 = boundary_edges.popitem()
            boundary.append(i_v0_init)
            boundary.append(i_v1)
            i_v0 = i_v1

            while True:
                try:
                    i_v1 = boundary_edges.pop(i_v0)
                    boundary.append(i_v1)
                    i_v0 = i_v1
                except KeyError:
                    if boundary[0] != boundary[-1]:
                        print('Boundary is not closed !!!')
                    else:
                        boundaries.append(boundary)
                    break
        except KeyError:
            break

    return {'v_v': v_v,
            'v_f': v_f,
            'f_f': f_f,
            'boundaries': boundaries}

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
