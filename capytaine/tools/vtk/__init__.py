#!/usr/bin/env python
# coding: utf-8
"""
"""


def compute_node_data(mesh, face_data):
    """Transform data defined at the center of the faces to data defined at the nodes of the mesh
    by a simple averaging of the values of the neighboring faces.

    Parameters
    ----------
    mesh: Mesh or CollectionOfMeshes
        the mesh on which the face face_data are defined
    face_data: numpy array of shape (mesh.nb_faces, ...)
        the data defined on the center of the faces of the mesh

    Returns
    -------
    node_data: numpy array of shape (mesh.nb_vertices, ...)
        the same data averaged on the nodes
    """

    import numpy as np
    from meshmagick.mesh import Mesh

    mesh = mesh if isinstance(mesh, Mesh) else mesh.merge()
    assert face_data.shape[0] == mesh.nb_faces

    # Initialize output array
    node_data_shape = (mesh.vertices.shape[0], ) + face_data.shape[1:]
    node_data = np.zeros(node_data_shape, dtype=np.complex)

    # Keep track of the number of faces near each vertex
    faces_near_nodes_shape = (mesh.vertices.shape[0], ) + (1, ) * len(face_data.shape[1:])
    faces_near_nodes = np.zeros(faces_near_nodes_shape, dtype=np.int8)

    for i, vertices in enumerate(mesh.faces):
        for vertex in vertices:
            faces_near_nodes[vertex] += 1
            node_data[vertex, ...] += face_data[i, ...]

    node_data /= faces_near_nodes
    return node_data

