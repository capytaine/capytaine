"""Importing mesh from meshio"""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
import logging

import numpy as np

from capytaine.tools.optional_imports import import_optional_dependency
from capytaine.meshes.meshes import Mesh

LOG = logging.getLogger(__name__)

def load_from_meshio(mesh, name=None):
    """Create a Mesh from a meshio mesh object."""
    meshio = import_optional_dependency("meshio")
    if not isinstance(mesh, meshio._mesh.Mesh):
        raise TypeError('mesh must be of type meshio._mesh.Mesh, received {:}'.format(type(mesh)))

    def all_faces_as_quads(cells):
        all_faces = []
        if 'quad' in cells:
            all_faces.append(cells['quad'])
        if 'triangle' in cells:
            num_triangles = len(mesh.cells_dict['triangle'])
            LOG.info("Stored {:} triangle faces as quadrilaterals".format(num_triangles))
            triangles_as_quads = np.empty((cells['triangle'].shape[0], 4), dtype=int)
            triangles_as_quads[:, :3] = cells['triangle'][:, :]
            triangles_as_quads[:, 3] = cells['triangle'][:, 2]  # Repeat one node to make a quad
            all_faces.append(triangles_as_quads)
        return np.concatenate(all_faces)

    if name is None:
        name = f'mesh_from_meshio_{next(Mesh._ids)}'

    mesh = Mesh(vertices=mesh.points, faces=all_faces_as_quads(mesh.cells_dict), name=name)
    mesh.heal_mesh()

    return mesh
