"""This module implements a tools to clip meshes against a plane.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np

from capytaine.meshes.geometry import Plane
from capytaine.meshes.meshes import Mesh

LOG = logging.getLogger(__name__)


def clip(source_mesh: Mesh, plane: Plane, vicinity_tol=1e-12, name=None):
    """Return a new mesh containing the source mesh clipped by the plane.

    Parameters
    ----------
    source_mesh : Mesh
        The mesh to be clipped.
    plane : Plane, optional
        The clipping plane.
    vicinity_tol : float, optional
        The absolute tolerance to consider en vertex is on the plane. Default is 1e-12.
    name: string, optional
        A name for the new clipped mesh.
    """
    vertices_data = _vertices_positions_wrt_plane(source_mesh, plane, vicinity_tol)

    nb_vertices_strictly_above_plane = np.count_nonzero(
        vertices_data['vertices_above_mask']
    )
    nb_vertices_below_or_on_plane = np.count_nonzero(
        vertices_data['vertices_below_mask'] | vertices_data['vertices_on_mask']
    )

    if nb_vertices_strictly_above_plane == source_mesh.nb_vertices:
        LOG.warning(f"Clipping {source_mesh.name} by {plane}: all vertices are removed.")
        clipped_mesh = Mesh(None, None)
        clipped_mesh._clipping_data = dict(faces_ids=[])

    elif nb_vertices_below_or_on_plane == source_mesh.nb_vertices:
        LOG.info(f"Clipping {source_mesh.name} by {plane}: no action.")
        clipped_mesh = source_mesh.copy()
        clipped_mesh._clipping_data = dict(faces_ids=list(range(source_mesh.nb_faces)))

    else:
        upper_mesh, crown_mesh, lower_mesh = _partition_mesh(vertices_data, source_mesh)

        if crown_mesh.nb_faces > 0:
            clipped_crown_mesh = _clip_crown(crown_mesh, plane)
            clipped_mesh = lower_mesh + clipped_crown_mesh
            clipped_mesh._clipping_data = {
                'faces_ids': np.concatenate((lower_mesh._clipping_data['faces_ids'],
                                             clipped_crown_mesh._clipping_data['faces_ids']))
            }
        else:
            clipped_mesh = lower_mesh

    if name is None:
        clipped_mesh.name = f'{source_mesh.name}_clipped'
    clipped_mesh.remove_unused_vertices()

    return clipped_mesh


def _vertices_positions_wrt_plane(source_mesh, plane, vicinity_tol):
    """Classifies vertices with respect to the clipping plane."""
    vertices_distances = plane.distance_to_point(source_mesh.vertices)
    vertices_data = {'vertices_distances': vertices_distances,
                     'vertices_above_mask': vertices_distances > vicinity_tol,
                     'vertices_on_mask': np.abs(vertices_distances) < vicinity_tol,
                     'vertices_below_mask': vertices_distances < -vicinity_tol
                     }
    return vertices_data


def _partition_mesh(vertices_data, source_mesh):
    """Partitions the mesh in 3 with respect to the plane:
    * upper_mesh: part entirely above the clipping plane
    * crown_mesh: part intersecting the clipping plane
    * lower_mesh: part entirely under the clipping plane
    """
    nb_vertices_above_per_face = vertices_data['vertices_above_mask'][source_mesh.faces].sum(axis=1)
    nb_vertices_below_per_face = vertices_data['vertices_below_mask'][source_mesh.faces].sum(axis=1)

    # Simple criteria ensuring that _faces are totally above or below the plane (4 _vertices at the same side)
    # Works for both triangles and quadrangles
    above_faces_mask = nb_vertices_above_per_face == 4
    below_faces_mask = nb_vertices_below_per_face == 4
    crown_faces_mask = np.logical_not(np.logical_or(above_faces_mask, below_faces_mask))

    faces_ids = {
        'upper_mesh': np.where(above_faces_mask)[0],
        'crown_mesh': np.where(crown_faces_mask)[0],
        'lower_mesh': np.where(below_faces_mask)[0]
    }

    partition = []
    for name in ['upper_mesh', 'crown_mesh', 'lower_mesh']:
        new_mesh, vertices_ids = source_mesh.extract_faces(
            id_faces_to_extract=faces_ids[name], return_index=True, name=name)
        new_mesh._clipping_data = {
            'faces_ids': faces_ids[name],
            'vertices_ids': vertices_ids,
            'vertices_distances': vertices_data['vertices_distances'][vertices_ids],
            'above_vertices_mask': vertices_data['vertices_above_mask'][vertices_ids],
            'on_vertices_mask': vertices_data['vertices_on_mask'][vertices_ids],
            'below_vertices_mask': vertices_data['vertices_below_mask'][vertices_ids],
        }
        partition.append(new_mesh)

    return partition


def _clip_crown(crown_mesh, plane):
    """Performs the clipping operation of the crown_mesh and determines the obtained boundaries.
    This is the heart method of the class.
    """
    vertices = crown_mesh.vertices

    above_vertices_mask = crown_mesh._clipping_data['above_vertices_mask']
    below_vertices_mask = crown_mesh._clipping_data['below_vertices_mask']

    on_vertices_mask = crown_mesh._clipping_data['on_vertices_mask']
    # TODO: Vertices pre-projection to be done here

    vertices_distances = crown_mesh._clipping_data['vertices_distances']

    # Init
    clipped_crown_mesh_faces = list()
    clipped_crown_mesh_relative_faces_ids = list()

    direct_boundary_edges = dict()
    inv_boundary_edges = dict()
    intersections_vertices = list()

    index_new_vertices = crown_mesh.nb_vertices

    for face_id in range(crown_mesh.nb_faces):

        face = crown_mesh.get_face(face_id)

        # # Determining the type of face clipping
        v_above_face = np.where(above_vertices_mask[face])[0]
        v_on_face = np.where(on_vertices_mask[face])[0]
        v_below_face = np.where(below_vertices_mask[face])[0]

        nb_above = len(v_above_face)
        nb_on = len(v_on_face)
        nb_below = len(v_below_face)

        face_type = str(nb_above) + str(nb_on) + str(nb_below)

        if face_type == '202':
            #    0*-----*3
            #     |     |
            # ----o-----o----
            #     |     |
            #    1*-----*2
            if v_above_face[1] == v_above_face[0] + 1:
                face = np.roll(face, -v_above_face[1])
            p0, p1, p2, p3 = vertices[face]
            ileft = plane.get_edge_intersection(p0, p1)
            iright = plane.get_edge_intersection(p2, p3)
            intersections_vertices += [ileft, iright]
            boundary_edge = [index_new_vertices, index_new_vertices + 1]
            clipped_crown_mesh_faces.append([index_new_vertices, face[1], face[2], index_new_vertices + 1])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 2

        elif face_type == '301':
            #      *2
            #     / \
            #    /   \
            #   /     \
            # 3*       *1
            #   \     /
            # ---o---o---
            #     \ /
            #      *0
            face = np.roll(face, -v_below_face[0])
            p0, p1, p3 = vertices[face[[0, 1, 3]]]
            ileft = plane.get_edge_intersection(p0, p3)
            iright = plane.get_edge_intersection(p0, p1)
            intersections_vertices += [ileft, iright]
            boundary_edge = [index_new_vertices, index_new_vertices + 1]
            clipped_crown_mesh_faces.append([index_new_vertices, face[0], index_new_vertices + 1, index_new_vertices])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 2

        elif face_type == '103':
            #      *0
            #     / \
            # ---o---o---
            #   /     \
            # 1* - - - *3
            #   \     /
            #    \   /
            #     \ /
            #      *2
            face = np.roll(face, -v_above_face[0])
            p0, p1, p3 = vertices[face[[0, 1, 3]]]
            ileft = plane.get_edge_intersection(p0, p1)
            iright = plane.get_edge_intersection(p0, p3)
            intersections_vertices += [ileft, iright]
            boundary_edge = [index_new_vertices, index_new_vertices + 1]
            clipped_crown_mesh_faces.append([index_new_vertices, face[1], face[3], index_new_vertices + 1])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            clipped_crown_mesh_faces.append([face[1], face[2], face[3], face[1]])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 2

        elif face_type == '102':
            #      *O
            #     / \
            # ---o---o---
            #   /     \
            # 1*-------*2
            face = np.roll(face, -v_above_face[0])
            p0, p1, p2 = vertices[face]
            ileft = plane.get_edge_intersection(p0, p1)
            iright = plane.get_edge_intersection(p0, p2)
            intersections_vertices += [ileft, iright]
            boundary_edge = [index_new_vertices, index_new_vertices + 1]
            clipped_crown_mesh_faces.append([index_new_vertices, face[1], face[2], index_new_vertices + 1])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 2

        elif face_type == '201':  # done
            #  2*-------*1
            #    \     /
            #  ---o---o---
            #      \ /
            #       *0
            face = np.roll(face, -v_below_face[0])
            p0, p1, p2 = vertices[face]
            ileft = plane.get_edge_intersection(p0, p2)
            iright = plane.get_edge_intersection(p0, p1)
            intersections_vertices += [ileft, iright]
            boundary_edge = [index_new_vertices, index_new_vertices + 1]
            clipped_crown_mesh_faces.append([index_new_vertices, face[0], index_new_vertices + 1, index_new_vertices])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 2

        elif face_type == '211':
            #        *3                   *1
            #       / \                  / \
            #      /   *2       or     2*   \
            #    0/   /                  \   \0
            # ---*---o---              ---o---*---
            #     \ /                      \ /
            #      *1                       *3
            #

            face = np.roll(face, -v_on_face[0])
            if vertices_distances[face[1]] < 0.:
                p1, p2 = vertices[face[[1, 2]]]
                iright = plane.get_edge_intersection(p1, p2)
                intersections_vertices.append(iright)
                boundary_edge = [face[0], index_new_vertices]
                clipped_crown_mesh_faces.append([face[0], face[1], index_new_vertices, face[0]])
            else:
                p2, p3 = vertices[face[[2, 3]]]
                ileft = plane.get_edge_intersection(p2, p3)
                intersections_vertices.append(ileft)
                boundary_edge = [index_new_vertices, face[0]]
                clipped_crown_mesh_faces.append([index_new_vertices, face[3], face[0], index_new_vertices])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 1

        elif face_type == '112':
            #       *3                     *1
            #      / \                    / \
            #  ---*---o---      or    ---o---*---
            #     0\   \                /   /0
            #       \   *2            2*   /
            #        \ /                \ /
            #         *1                 *3
            face = np.roll(face, -v_on_face[0])
            if vertices_distances[face[1]] < 0.:
                p2, p3 = vertices[face[[2, 3]]]
                iright = plane.get_edge_intersection(p2, p3)
                intersections_vertices.append(iright)
                boundary_edge = [face[0], index_new_vertices]
                clipped_crown_mesh_faces.append([face[0], face[1], face[2], index_new_vertices])
            else:
                p1, p2 = vertices[face[[1, 2]]]
                ileft = plane.get_edge_intersection(p1, p2)
                intersections_vertices.append(ileft)
                boundary_edge = [index_new_vertices, face[0]]
                clipped_crown_mesh_faces.append([index_new_vertices, face[2], face[3], face[0]])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 1

        elif face_type == '013':
            # -----*-----
            #     / \
            #    /   \
            #   *     *
            #    \   /
            #     \ /
            #      *
            boundary_edge = None
            clipped_crown_mesh_faces.append(list(face))
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '210' or face_type == '310':
            #   *-------*               *
            #    \ 210 /               / \ 310
            #     \   /               *   *
            #      \ /                 \ /
            #   ----*----           ----*----
            boundary_edge = None

        elif face_type == '111':
            #        *2              *1
            #       /|               |\
            #      / |               | \
            #  ---*--o---    or   ---o--*---
            #     0\ |               | /0
            #       \|               |/
            #        *1              *2
            face = np.roll(face, -v_on_face[0])
            p1, p2 = vertices[face[[1, 2]]]
            if vertices_distances[face[1]] < 0.:
                iright = plane.get_edge_intersection(p1, p2)
                intersections_vertices.append(iright)
                boundary_edge = [face[0], index_new_vertices]
                clipped_crown_mesh_faces.append([face[0], face[1], index_new_vertices, face[0]])
            else:
                ileft = plane.get_edge_intersection(p1, p2)
                intersections_vertices.append(ileft)
                boundary_edge = [index_new_vertices, face[0]]
                clipped_crown_mesh_faces.append([index_new_vertices, face[2], face[0], index_new_vertices])
            clipped_crown_mesh_relative_faces_ids.append(face_id)
            index_new_vertices += 1

        elif face_type == '120':
            #         *O
            #        / \
            #       /   \
            #     1/     \2
            # ----*-------*----
            # face = np.roll(face, -v_above_face[0])
            # boundary_edge = [face[1], face[2]]
            # FIXME: quick fix here : robust ?
            boundary_edge = None

        elif face_type == '021':
            #  ----*-------*----
            #      2\     /1
            #        \   /
            #         \ /
            #          *0
            face = np.roll(face, -v_below_face[0])
            boundary_edge = [face[2], face[1]]
            face = list(face)
            face.append(face[0])
            clipped_crown_mesh_faces.append(face)
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '022':
            # ----*-----*----
            #    0|     |3
            #     |     |
            #    1*-----*2
            if v_on_face[1] == v_on_face[0] + 1:
                face = np.roll(face, -v_on_face[1])
            boundary_edge = [face[0], face[3]]
            clipped_crown_mesh_faces.append(list(face))
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '012':
            #   ------*------
            #        / \
            #       /   \
            #      /     \
            #     *-------*
            boundary_edge = None
            face = list(face)
            face.append(face[0])
            clipped_crown_mesh_faces.append(face)
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '220':
            #    0*-----*3
            #     |     |
            #    1|     |2
            # ----*-----*----

            # if v_above_face[1] == v_above_face[0] + 1:
            #     face = np.roll(face, -v_above_face[1])
            # boundary_edge = [face[1], face[2]]
            # FIXME: quick fix here : robust ?
            boundary_edge = None

        elif face_type == '121':
            #       *0
            #      / \
            #     /   \
            # ---*-----*---
            #    1\   /3
            #      \ /
            #       *2
            face = np.roll(face, -v_above_face[0])
            boundary_edge = [face[1], face[3]]
            clipped_crown_mesh_faces.append([face[1], face[2], face[3], face[1]])
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '300' or face_type == '400':
            #       *               *-----*
            #      / \              |     |
            #     /300\      or     | 400 |
            #    *-----*            *-----*
            # ____________       ______________
            boundary_edge = None

        elif face_type == '003':
            #  -----------
            #       *
            #      / \
            #     /   \
            #    *-----*
            boundary_edge = None
            face = list(face)
            face.append(face[0])
            clipped_crown_mesh_faces.append(face)
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '004':
            #  ---------------
            #      *-----*
            #      |     |
            #      |     |
            #      *-----*
            boundary_edge = None
            clipped_crown_mesh_faces.append(list(face))
            clipped_crown_mesh_relative_faces_ids.append(face_id)

        elif face_type == '030' or face_type == '040':
            # Face is totally on the plane --> rare case...
            boundary_edge = None

        else:
            raise Exception("Face %u clipping case %s not known." % (face_id, face_type))

        # Building boundary connectivity
        if boundary_edge is not None:
            direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
            inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]

    if len(intersections_vertices) > 0:
        vertices = np.concatenate((vertices, intersections_vertices))

    clipped_crown_mesh = Mesh(vertices, clipped_crown_mesh_faces)
    clipped_crown_mesh._clipping_data = {
        'faces_ids': crown_mesh._clipping_data['faces_ids'][clipped_crown_mesh_relative_faces_ids]
    }

    return clipped_crown_mesh
