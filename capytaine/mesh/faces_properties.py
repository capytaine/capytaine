#!/usr/bin/env python
# coding: utf-8
"""Helper functions to compute the properties of the faces in the mesh."""

import numpy as np


def compute_faces_properties(mesh):
    """Updates the faces properties of the mesh"""

    # faces_areas, faces_normals, faces_centers = mm.get_all_faces_properties(mesh._vertices, mesh._faces)
    nf = mesh.nb_faces

    # triangle_mask = _faces[:, 0] == _faces[:, -1]
    # nb_triangles = np.sum(triangle_mask)
    # quads_mask = np.invert(triangle_mask)
    # nb_quads = nf - nb_triangles

    faces_areas = np.zeros(nf, dtype=np.float)
    faces_normals = np.zeros((nf, 3), dtype=np.float)
    faces_centers = np.zeros((nf, 3), dtype=np.float)

    # Collectively dealing with triangles
    # triangles = _faces[triangle_mask]
    triangles_id = mesh.triangles_ids
    triangles = mesh._faces[triangles_id]

    triangles_normals = np.cross(mesh._vertices[triangles[:, 1]] - mesh._vertices[triangles[:, 0]],
                                 mesh._vertices[triangles[:, 2]] - mesh._vertices[triangles[:, 0]])
    triangles_areas = np.linalg.norm(triangles_normals, axis=1)
    faces_normals[triangles_id] = triangles_normals / np.array(([triangles_areas, ] * 3)).T
    faces_areas[triangles_id] = triangles_areas / 2.
    faces_centers[triangles_id] = np.sum(mesh._vertices[triangles[:, :3]], axis=1) / 3.

    # Collectively dealing with quads
    quads_id = mesh.quadrangles_ids
    quads = mesh._faces[quads_id]
    # quads = _faces[quads_mask]

    quads_normals = np.cross(mesh._vertices[quads[:, 2]] - mesh._vertices[quads[:, 0]],
                             mesh._vertices[quads[:, 3]] - mesh._vertices[quads[:, 1]])
    faces_normals[quads_id] = quads_normals / np.array(([np.linalg.norm(quads_normals, axis=1), ] * 3)).T

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

