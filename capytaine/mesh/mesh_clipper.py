#!/usr/bin/env python
# coding: utf-8
"""This module holds a tools to clip meshes against a plane.
Based on Meshmagick by François Rongère (EC Nantes).
"""

import logging

import numpy as np

from capytaine.mesh.mesh import Mesh

LOG = logging.getLogger(__name__)


def clip(source_mesh, plane, vicinity_tol=1e-3, return_all_data=False):
    """
    Parameters
    ----------
    source_mesh : Mesh
        The mesh to be clipped.
    plane : Plane, optional
        The clipping plane. By default, the plane is the Oxy plane.
    vicinity_tol : float, optional
        The absolute tolerance to consider en vertex is on the plane. Default is 1e-3.
    return_all_data: bool, optional
        If False, return only the clipped mesh, else returned detailed data
    """
    data = {}
    data.update(_vertices_positions_wrt_plane(source_mesh, plane, vicinity_tol))

    if source_mesh.nb_vertices == data['nb_vertices_above_or_on_mask']:
        LOG.warning(f"Clipping {source_mesh.name} by {plane}: "
                    "all vertices are removed.")
        data['clipped_mesh'] = Mesh(None, None)
        data['clipped_mesh_faces_ids'] = []

    elif source_mesh.nb_vertices == data['nb_vertices_below_or_on_mask']:
        LOG.info(f"Clipping {source_mesh.name} by {plane}: "
                 "no action.")
        data['clipped_mesh'] = source_mesh
        data['clipped_mesh_faces_ids'] = list(range(source_mesh.nb_faces))

    else:
        data.update(_partition_mesh(data, source_mesh))

        if data['crown_mesh'].nb_faces > 0:
            data.update(_clip_crown_by_plane(data, source_mesh, plane))
            clipped_mesh = data['lower_mesh'] + data['clipped_crown_mesh']
            data['clipped_mesh_faces_ids'] = list(data['below_faces_ids']) + list(data['clipped_crown_mesh_faces_ids'])
        else:
            clipped_mesh = data['lower_mesh']
            data['clipped_mesh_faces_ids'] = list(data['below_faces_ids'])

        data['clipped_mesh'] = clipped_mesh

    clipped_mesh_name = '_'.join((source_mesh.name, 'clipped'))
    data['clipped_mesh'].name = clipped_mesh_name
    data['clipped_mesh'].remove_unused_vertices()

    if return_all_data:
        return data
    else:
        return data['clipped_mesh']


def _vertices_positions_wrt_plane(source_mesh, plane, vicinity_tol):
    """
    Classifies vertices with respect to the clipping plane
    """
    vertices_distances = plane.distance_to_point(source_mesh.vertices)

    vertices_positions = {'vertices_distances': vertices_distances,
                          'vertices_above_mask': vertices_distances > vicinity_tol,
                          'vertices_on_mask': np.fabs(vertices_distances) < vicinity_tol,
                          'vertices_below_mask': vertices_distances < -vicinity_tol
                          }
    nb_vertices = {
        'nb_vertices_above_or_on_mask': np.count_nonzero(
            vertices_positions['vertices_above_mask']| vertices_positions['vertices_on_mask']
        ),
        'nb_vertices_below_or_on_mask': np.count_nonzero(
            vertices_positions['vertices_below_mask'] | vertices_positions['vertices_on_mask']
        ),
    }
    vertices_positions.update(nb_vertices)
    return vertices_positions


def _partition_mesh(data, source_mesh):
    """Partitions the mesh in 3 with respect to the plane

    * upper_mesh: part entirely above the clipping plane
    * crown_mesh: part intersecting the clipping plane
    * lower_mesh: part entirely under the clipping plane
    """

    vertices_distances = data['vertices_distances']
    vertices_above_mask = data['vertices_above_mask']
    vertices_on_mask = data['vertices_on_mask']
    vertices_below_mask = data['vertices_below_mask']

    source_mesh_faces = source_mesh.faces

    nb_vertices_above = vertices_above_mask[source_mesh_faces].sum(axis=1)
    nb_vertices_below = vertices_below_mask[source_mesh_faces].sum(axis=1)

    # Simple criteria ensuring that _faces are totally above or below the plane (4 _vertices at the same side)
    # Works for both triangles and quadrangles
    above_faces_mask = nb_vertices_above == 4
    below_faces_mask = nb_vertices_below == 4
    crown_faces_mask = np.logical_and(np.logical_not(above_faces_mask), np.logical_not(below_faces_mask))

    above_faces_ids = np.where(above_faces_mask)[0]
    below_faces_ids = np.where(below_faces_mask)[0]
    crown_faces_ids = np.where(crown_faces_mask)[0]

    partition = dict()

    def generate_mesh(faces_ids, key):
        new_mesh, ids = source_mesh.extract_faces(faces_ids, return_index=True)
        new_mesh.name = key
        partition['_'.join((key, 'vertices_ids'))] = ids
        partition['_'.join((key, 'vertices_distances'))] = vertices_distances[ids]
        partition['_'.join((key, 'above_vertices_mask'))] = vertices_above_mask[ids]
        partition['_'.join((key, 'on_vertices_mask'))] = vertices_on_mask[ids]
        partition['_'.join((key, 'below_vertices_mask'))] = vertices_below_mask[ids]
        partition[key] = new_mesh

    # Generating partition meshes
    generate_mesh(above_faces_ids, 'upper_mesh')
    generate_mesh(crown_faces_ids, 'crown_mesh')
    generate_mesh(below_faces_ids, 'lower_mesh')
    partition['above_faces_ids'] = above_faces_ids
    partition['crown_faces_ids'] = crown_faces_ids
    partition['below_faces_ids'] = below_faces_ids

    return partition


def _clip_crown_by_plane(data, source_mesh, plane):
    """Performs the clipping operation of the crown_mesh and determines the obtained boundaries.

    This is the heart method of the class.
    """

    crown_mesh = data['crown_mesh']
    vertices = crown_mesh.vertices
    vertices_on_mask = data['crown_mesh_on_vertices_mask']

    # TODO: Vertices pre-projection to be done here !!!
    # vertices_on = partition['vertices_on']
    # _vertices[vertices_on] = plane.orthogonal_projection_on_plane(_vertices[vertices_on])
    # pos[vertices_on] = 0.

    vertices_above_mask = data['crown_mesh_above_vertices_mask']
    vertices_below_mask = data['crown_mesh_below_vertices_mask']

    vertices_distances = data['crown_mesh_vertices_distances']

    # Init
    crown_faces = list()
    clipped_crown_relative_faces_ids = list()
    direct_boundary_edges = dict()
    inv_boundary_edges = dict()
    intersections = list()

    index = crown_mesh.nb_vertices

    for face_id in range(crown_mesh.nb_faces):

        face = crown_mesh.get_face(face_id)

        # # Determining the type of face clipping
        v_above_face = np.where(vertices_above_mask[face])[0]
        v_on_face = np.where(vertices_on_mask[face])[0]
        v_below_face = np.where(vertices_below_mask[face])[0]

        nb_above = len(v_above_face)
        nb_on = len(v_on_face)
        nb_below = len(v_below_face)

        face_type = str(nb_above) + str(nb_on) + str(nb_below)

        if face_type == '202':  # Done
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
            intersections += [ileft, iright]
            boundary_edge = [index, index + 1]
            crown_faces.append([index, face[1], face[2], index + 1])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 2

        elif face_type == '301':  # Done
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
            intersections += [ileft, iright]
            boundary_edge = [index, index + 1]
            crown_faces.append([index, face[0], index + 1, index])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 2

        elif face_type == '103':  # Done
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
            intersections += [ileft, iright]
            boundary_edge = [index, index + 1]
            crown_faces.append([index, face[1], face[3], index + 1])
            clipped_crown_relative_faces_ids.append(face_id)
            crown_faces.append([face[1], face[2], face[3], face[1]])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 2

        elif face_type == '102':  # Done
            #      *O
            #     / \
            # ---o---o---
            #   /     \
            # 1*-------*2
            face = np.roll(face, -v_above_face[0])
            p0, p1, p2 = vertices[face]
            ileft = plane.get_edge_intersection(p0, p1)
            iright = plane.get_edge_intersection(p0, p2)
            intersections += [ileft, iright]
            boundary_edge = [index, index + 1]
            crown_faces.append([index, face[1], face[2], index + 1])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 2

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
            intersections += [ileft, iright]
            boundary_edge = [index, index + 1]
            crown_faces.append([index, face[0], index + 1, index])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 2

        elif face_type == '211':  # Done
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
                intersections.append(iright)
                boundary_edge = [face[0], index]
                crown_faces.append([face[0], face[1], index, face[0]])
            else:
                p2, p3 = vertices[face[[2, 3]]]
                ileft = plane.get_edge_intersection(p2, p3)
                intersections.append(ileft)
                boundary_edge = [index, face[0]]
                crown_faces.append([index, face[3], face[0], index])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 1

        elif face_type == '112':  # Done
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
                intersections.append(iright)
                boundary_edge = [face[0], index]
                crown_faces.append([face[0], face[1], face[2], index])
            else:
                p1, p2 = vertices[face[[1, 2]]]
                ileft = plane.get_edge_intersection(p1, p2)
                intersections.append(ileft)
                boundary_edge = [index, face[0]]
                crown_faces.append([index, face[2], face[3], face[0]])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 1

        elif face_type == '013':  # Done
            # -----*-----
            #     / \
            #    /   \
            #   *     *
            #    \   /
            #     \ /
            #      *
            boundary_edge = None
            crown_faces.append(list(face))
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '210' or face_type == '310':  # Done
            #   *-------*               *
            #    \ 210 /               / \ 310
            #     \   /               *   *
            #      \ /                 \ /
            #   ----*----           ----*----
            boundary_edge = None

        elif face_type == '111':  # Done
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
                intersections.append(iright)
                boundary_edge = [face[0], index]
                crown_faces.append([face[0], face[1], index, face[0]])
            else:
                ileft = plane.get_edge_intersection(p1, p2)
                intersections.append(ileft)
                boundary_edge = [index, face[0]]
                crown_faces.append([index, face[2], face[0], index])
            clipped_crown_relative_faces_ids.append(face_id)
            index += 1

        elif face_type == '120':  # Done
            #         *O
            #        / \
            #       /   \
            #     1/     \2
            # ----*-------*----
            # face = np.roll(face, -v_above_face[0])
            # boundary_edge = [face[1], face[2]]
            # FIXME: quick fix here : robust ?
            boundary_edge = None

        elif face_type == '021':  # Done
            #  ----*-------*----
            #      2\     /1
            #        \   /
            #         \ /
            #          *0
            face = np.roll(face, -v_below_face[0])
            boundary_edge = [face[2], face[1]]
            face = list(face)
            face.append(face[0])
            crown_faces.append(face)
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '022':
            # ----*-----*----
            #    0|     |3
            #     |     |
            #    1*-----*2
            if v_on_face[1] == v_on_face[0] + 1:
                face = np.roll(face, -v_on_face[1])
            boundary_edge = [face[0], face[3]]
            crown_faces.append(list(face))
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '012':  # Done
            #   ------*------
            #        / \
            #       /   \
            #      /     \
            #     *-------*
            boundary_edge = None
            face = list(face)
            face.append(face[0])
            crown_faces.append(face)
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '220':  # Done
            #    0*-----*3
            #     |     |
            #    1|     |2
            # ----*-----*----

            # if v_above_face[1] == v_above_face[0] + 1:
            #     face = np.roll(face, -v_above_face[1])
            # boundary_edge = [face[1], face[2]]
            # FIXME: quick fix here : robust ?
            boundary_edge = None

        elif face_type == '121':  # Done
            #       *0
            #      / \
            #     /   \
            # ---*-----*---
            #    1\   /3
            #      \ /
            #       *2
            face = np.roll(face, -v_above_face[0])
            boundary_edge = [face[1], face[3]]
            crown_faces.append([face[1], face[2], face[3], face[1]])
            clipped_crown_relative_faces_ids.append(face_id)

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
            crown_faces.append(face)
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '004':
            #  ---------------
            #      *-----*
            #      |     |
            #      |     |
            #      *-----*
            boundary_edge = None
            crown_faces.append(list(face))
            clipped_crown_relative_faces_ids.append(face_id)

        elif face_type == '030' or face_type == '040':
            # Face is totally on the plane --> rare case...
            boundary_edge = None

        else:
            try:
                from capytaine.io import mesh_writers
                mesh_writers.write_VTP('full_debug.vtp', source_mesh.vertices, source_mesh.faces)
                # meshio.write_VTP('clipped_crown_debug.vtp', clipped_crown_mesh.vertices, clipped_crown_mesh.faces)
                mesh_writers.write_VTP('crown_debug.vtp', crown_mesh.vertices, crown_mesh.faces)
            except:
                pass
            raise Exception("Face %u clipping case %s not known." % (face_id, face_type))

        # Building boundary connectivity
        if boundary_edge is not None:
            direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
            inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]

    if len(intersections) > 0:
        vertices = np.concatenate((vertices, intersections))

    clipped_crown_mesh = Mesh(vertices, crown_faces)

    return {'clipped_crown_mesh': clipped_crown_mesh,
            'clipped_crown_mesh_faces_ids': data['crown_faces_ids'][clipped_crown_relative_faces_ids]}

