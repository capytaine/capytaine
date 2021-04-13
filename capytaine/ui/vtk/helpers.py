#!/usr/bin/env python
# coding: utf-8
"""Tools for 3D displays with VTK."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

from typing import Union

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.tools.optional_imports import import_optional_dependency

vtk = import_optional_dependency("vtk")

def compute_vtk_polydata(mesh: Union[Mesh, CollectionOfMeshes]):
    """Transform a mesh into vtkPolydata."""

    # Create a vtkPoints object and store the points in it
    points = vtk.vtkPoints()
    for point in mesh.vertices:
        points.InsertNextPoint(point)

    # Create a vtkCellArray to store faces
    faces = vtk.vtkCellArray()
    for face_ids in mesh.faces:
        if face_ids[0] == face_ids[-1]:
            # Triangle
            curface = face_ids[:3]
            vtk_face = vtk.vtkTriangle()
        else:
            # Quadrangle
            curface = face_ids[:4]
            vtk_face = vtk.vtkQuad()

        for idx, id in enumerate(curface):
            vtk_face.GetPointIds().SetId(idx, id)

        faces.InsertNextCell(vtk_face)

    vtk_polydata = vtk.vtkPolyData()
    vtk_polydata.SetPoints(points)
    vtk_polydata.SetPolys(faces)

    return vtk_polydata


def compute_node_data(mesh: Union[Mesh, CollectionOfMeshes],
                      face_data):
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

    mesh = mesh.merged()
    assert face_data.shape[0] == mesh.nb_faces

    # Initialize output array
    node_data_shape = (mesh.vertices.shape[0], ) + face_data.shape[1:]
    node_data = np.zeros(node_data_shape, dtype=complex)

    # Keep track of the number of faces near each vertex
    faces_near_nodes_shape = (mesh.vertices.shape[0], ) + (1, ) * len(face_data.shape[1:])
    nb_faces_near_nodes = np.zeros(faces_near_nodes_shape, dtype=np.int8)

    for i, vertices in enumerate(mesh.faces):
        for vertex in vertices:
            nb_faces_near_nodes[vertex] += 1
            node_data[vertex, ...] += face_data[i, ...]

    node_data /= nb_faces_near_nodes
    return node_data

