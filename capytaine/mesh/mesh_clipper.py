#!/usr/bin/env python
# coding: utf-8

"""This module holds a tools to clip meshes against a plane.

Based on Meshmagick by Francois Rongere (EC Nantes).
"""

from capytaine.mesh.mesh import *
from capytaine.tools.geometry import *


class MeshClipper:
    """A class to perform mesh clipping operations.

    Parameters
    ----------
    source_mesh : Mesh
        The mesh to be clipped.
    plane : Plane, optional
        The clipping plane. By default, the plane is the Oxy plane.
    vicinity_tol : float, optional
        The absolute tolerance to consider en vertex is on the plane. Default is 1e-3.
    assert_closed_boundaries : bool, optional
        False by default. When True, the mesh clipper will raise an exception if intersections with the clipping
        plane are not closed. It may be caused by a non-watertight mesh.
    verbose : bool, optional
        False by default. If True, some messages on operations that are handled are printed.
    """
    def __init__(self, source_mesh, plane=Plane(), vicinity_tol=1e-3, assert_closed_boundaries=False, verbose=False):
        self._source_mesh = source_mesh
        self._plane = plane

        self._vicinity_tol = vicinity_tol

        self._assert_closed_boundaries = assert_closed_boundaries
        self._verbose = verbose

        self.__internals__ = dict()

        self._update()

    @property
    def verbose(self):
        """Get the current verbosity"""

        return self._verbose

    def verbose_on(self):
        """Switches ON the verbosity of the clipper."""

        self._verbose = True

    def verbose_off(self):
        """Switches OFF the verbosity of the clipper."""

        self._verbose = False

    @property
    def assert_closed_boundaries(self):
        """Do we assert the boundaries have to be closed"""

        return self._assert_closed_boundaries

    def assert_closed_boundaries_on(self):
        """Switches ON the flag for closed boundary assertion while clipping."""

        self._assert_closed_boundaries = True

    def assert_closed_boundaries_off(self):
        """Switches OFF the flag for closed boundary assertion while clipping."""

        self._assert_closed_boundaries = False
        return

    @property
    def vicinity_tol(self):
        """Vicinity tolerance.

        It tells if a point is close enough to the plane to consider it lies on the plane
        """

        return self._vicinity_tol

    @vicinity_tol.setter
    def vicinity_tol(self, value):
        """Set the vicinity tolerance that tells that a vertex is located on the plane."""

        self.__internals__.clear()
        self._vicinity_tol = float(value)
        self._update()

    @property
    def source_mesh(self):
        """The mesh we work with"""

        return self._source_mesh

    @property
    def plane(self):
        """The clipping plane"""

        return self._plane

    @plane.setter
    def plane(self, value):
        """Changes the clipping plane."""

        self.__internals__.clear()
        self._plane = value
        self._update()

    def _update(self):
        """Updates the clipper"""

        self._vertices_positions_wrt_plane()
        self._partition_mesh()
        self._clip()

    def _vertices_positions_wrt_plane(self):
        """
        Classifies vertices with respect to the clipping plane
        """

        vertices_distances = self._plane.get_point_dist_wrt_plane(self._source_mesh.vertices)

        vertices_positions = {'vertices_distances': vertices_distances,
                              'vertices_above_mask': vertices_distances > self._vicinity_tol,
                              'vertices_on_mask': np.fabs(vertices_distances) < self._vicinity_tol,
                              'vertices_below_mask': vertices_distances < -self._vicinity_tol
                              }
        self.__internals__.update(vertices_positions)

    def _partition_mesh(self):
        """Partitions the mesh in 3 with respect to the plane

        * upper_mesh: part entirely above the clipping plane
        * crown_mesh: part intersecting the clipping plane
        * lower_mesh: part entirely under the clipping plane
        """

        vertices_distances = self.__internals__['vertices_distances']
        vertices_above_mask = self.__internals__['vertices_above_mask']
        vertices_on_mask = self.__internals__['vertices_on_mask']
        vertices_below_mask = self.__internals__['vertices_below_mask']

        source_mesh_faces = self.source_mesh.faces

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
            new_mesh, ids = self._source_mesh.extract_faces(faces_ids, return_index=True)
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

        self.__internals__.update(partition)

    @property
    def lower_mesh(self):
        """A new mesh composed of the faces that entirely lie under the clipping plane.

        Returns
        -------
        Mesh
        """

        return self.__internals__['lower_mesh']

    @property
    def crown_mesh(self):
        """A new mesh only having the faces that cut the plane

        Returns
        -------
        Mesh
        """

        return self.__internals__['crown_mesh']

    @property
    def upper_mesh(self):
        """A new mesh only having the faces lying entirely up the plane

        Returns
        -------
        Mesh
        """

        return self.__internals__['upper_mesh']

    @property
    def clipped_crown_mesh(self):
        """A new mesh that is obtained by clipping the crown_mesh.

        Returns
        -------
        Mesh
        """

        return self.__internals__['clipped_crown_mesh']

    @property
    def clipped_mesh(self):
        """The resulting clipped mesh"""

        return self.__internals__['clipped_mesh']

    @property
    def closed_polygons(self):
        """Returns the list of closed boundary polygons obtained after clipping.

        This is a list of lists. The enclosed lists are ordered IDs list that form a closed polygon, described in the
        counter-clockwise order with respect to the mesh (oriented following the clipping plane's normal).

        Returns
        -------
        list

        Warnings
        --------

        * The first vertex is repeated at the end of the list. By definition, these polygons are lying on the clipping
          plane.
        * Vertices IDs are corresponding to the IDs of the clipped_crown_mesh, not those of the clipped_mesh.
        """

        return self.__internals__['closed_polygons']

    @property
    def closed_polygons_vertices(self):
        """Returns the list of closed boundary polygons obtained after clipping.

        This is a list of lists. The enclosed lists are ordered vertices coordinates of the closed polygons. By
        definition, these polygons are lying on the clipping plane.

        Returns
        -------
        list
        """

        polygons = self.__internals__['closed_polygons']
        closed_polygons_vertices = []
        # TODO: voir si on ne peut pas directement indicer par polygons sans boucle for
        for polygon in polygons:
            closed_polygons_vertices.append(self.clipped_crown_mesh.vertices[polygon])
        return closed_polygons_vertices

    @property
    def nb_closed_polygons(self):
        """The number of closed polygons obtained after clipping

        Returns
        -------
        int
        """

        return len(self.__internals__['closed_polygons'])

    @property
    def open_lines(self):
        """Returns a list of open boundary lines obtained after clipping.

        This is a list of lists. The enclosed lists are ordered vertices IDs.

        Returns
        -------
        list


        Warning
        -------

        * The vertices IDs correspond to the IDs of clipped_crown_mesh, not clipped_mesh.
        """

        return self.__internals__['open_lines']

    @property
    def open_lines_vertices(self):
        """Returns a list of open boundary lines obtained after clipping.

        This is a list of lists. The enclosed lists are ordered vertices IDs.

        Returns
        -------
        list
        """

        lines = self.__internals__['open_lines']
        lines_vertices = []
        # TODO: voir si on ne peut pas directement indicer par polygons sans boucle for
        for line in lines:
            lines_vertices.append(self.clipped_crown_mesh.vertices[line])
        return lines_vertices

    @property
    def nb_open_lines(self):
        """The number of open lines obtained after clipping.

        Returns
        -------
        int
        """

        return len(self.__internals__['open_lines'])

    def _clip_crown_by_plane(self):
        """Performs the clipping operation of the crown_mesh and determines the obtained boundaries.

        This is the heart method of the class.
        """

        crown_mesh = self.crown_mesh
        vertices = crown_mesh.vertices
        vertices_on_mask = self.__internals__['crown_mesh_on_vertices_mask']

        # TODO: Vertices pre-projection to be done here !!!
        # vertices_on = partition['vertices_on']
        # _vertices[vertices_on] = plane.orthogonal_projection_on_plane(_vertices[vertices_on])
        # pos[vertices_on] = 0.

        vertices_above_mask = self.__internals__['crown_mesh_above_vertices_mask']
        vertices_below_mask = self.__internals__['crown_mesh_below_vertices_mask']

        vertices_distances = self.__internals__['crown_mesh_vertices_distances']

        # Init
        crown_faces = list()
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
                ileft = self._plane.get_edge_intersection(p0, p1)
                iright = self._plane.get_edge_intersection(p2, p3)
                intersections += [ileft, iright]
                boundary_edge = [index, index + 1]
                crown_faces.append([index, face[1], face[2], index + 1])
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
                ileft = self._plane.get_edge_intersection(p0, p3)
                iright = self._plane.get_edge_intersection(p0, p1)
                intersections += [ileft, iright]
                boundary_edge = [index, index + 1]
                crown_faces.append([index, face[0], index + 1, index])
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
                ileft = self._plane.get_edge_intersection(p0, p1)
                iright = self._plane.get_edge_intersection(p0, p3)
                intersections += [ileft, iright]
                boundary_edge = [index, index + 1]
                crown_faces.append([index, face[1], face[3], index + 1])
                crown_faces.append([face[1], face[2], face[3], face[1]])
                index += 2

            elif face_type == '102':  # Done
                #      *O
                #     / \
                # ---o---o---
                #   /     \
                # 1*-------*2
                face = np.roll(face, -v_above_face[0])
                p0, p1, p2 = vertices[face]
                ileft = self._plane.get_edge_intersection(p0, p1)
                iright = self._plane.get_edge_intersection(p0, p2)
                intersections += [ileft, iright]
                boundary_edge = [index, index + 1]
                crown_faces.append([index, face[1], face[2], index + 1])
                index += 2

            elif face_type == '201':  # done
                #  2*-------*1
                #    \     /
                #  ---o---o---
                #      \ /
                #       *0
                face = np.roll(face, -v_below_face[0])
                p0, p1, p2 = vertices[face]
                ileft = self._plane.get_edge_intersection(p0, p2)
                iright = self._plane.get_edge_intersection(p0, p1)
                intersections += [ileft, iright]
                boundary_edge = [index, index + 1]
                crown_faces.append([index, face[0], index + 1, index])
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
                    iright = self._plane.get_edge_intersection(p1, p2)
                    intersections.append(iright)
                    boundary_edge = [face[0], index]
                    crown_faces.append([face[0], face[1], index, face[0]])
                else:
                    p2, p3 = vertices[face[[2, 3]]]
                    ileft = self._plane.get_edge_intersection(p2, p3)
                    intersections.append(ileft)
                    boundary_edge = [index, face[0]]
                    crown_faces.append([index, face[3], face[0], index])
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
                    iright = self._plane.get_edge_intersection(p2, p3)
                    intersections.append(iright)
                    boundary_edge = [face[0], index]
                    crown_faces.append([face[0], face[1], face[2], index])
                else:
                    p1, p2 = vertices[face[[1, 2]]]
                    ileft = self._plane.get_edge_intersection(p1, p2)
                    intersections.append(ileft)
                    boundary_edge = [index, face[0]]
                    crown_faces.append([index, face[2], face[3], face[0]])
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
                    iright = self._plane.get_edge_intersection(p1, p2)
                    intersections.append(iright)
                    boundary_edge = [face[0], index]
                    crown_faces.append([face[0], face[1], index, face[0]])
                else:
                    ileft = self._plane.get_edge_intersection(p1, p2)
                    intersections.append(ileft)
                    boundary_edge = [index, face[0]]
                    crown_faces.append([index, face[2], face[0], index])
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

            elif face_type == '022':
                # ----*-----*----
                #    0|     |3
                #     |     |
                #    1*-----*2
                if v_on_face[1] == v_on_face[0] + 1:
                    face = np.roll(face, -v_on_face[1])
                boundary_edge = [face[0], face[3]]
                crown_faces.append(list(face))

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

            elif face_type == '004':
                #  ---------------
                #      *-----*
                #      |     |
                #      |     |
                #      *-----*
                boundary_edge = None
                crown_faces.append(list(face))

            elif face_type == '030' or face_type == '040':
                # Face is totally on the plane --> rare case...
                boundary_edge = None

            else:
                try:
                    from . import mmio
                    mmio.write_VTP('full_debug.vtp', self.source_mesh.vertices, self.source_mesh.faces)
                    # mmio.write_VTP('clipped_crown_debug.vtp', clipped_crown_mesh.vertices, clipped_crown_mesh.faces)
                    mmio.write_VTP('crown_debug.vtp', crown_mesh.vertices, crown_mesh.faces)
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

        # TODO: faire un merge uniquement sur la liste instersections et non sur tout le maillage clipped_crown
        # FIXME: potentiellement, un bug a ete introduit ici !!! --> l'update n'est plus bon sur les dictionnaires...
        new_id = clipped_crown_mesh.merge_duplicates(atol=1e-5)  # Warning: choosing a lower value

        # Updating dictionaries
        direct_boundary_edges = dict(
            list(zip(new_id[list(direct_boundary_edges.keys())], new_id[list(direct_boundary_edges.values())])))
        inv_boundary_edges = dict(list(zip(new_id[list(inv_boundary_edges.keys())], new_id[list(inv_boundary_edges.values())])))

        # Ordering boundary edges in continuous lines
        closed_polygons = list()
        open_lines = list()
        while True:
            try:
                line = list()
                v0_init, v1 = direct_boundary_edges.popitem()
                line.append(v0_init)
                line.append(v1)
                v0 = v1

                while True:
                    try:
                        v1 = direct_boundary_edges.pop(v0)
                        line.append(v1)
                        v0 = v1
                    except KeyError:
                        if line[0] != line[-1]:
                            # Trying to find an other queue
                            queue = list()
                            v0 = v0_init
                            while True:
                                try:
                                    v1 = inv_boundary_edges[v0]
                                    direct_boundary_edges.pop(v1)
                                    queue.append(v1)
                                    v0 = v1
                                except:
                                    queue.reverse()
                                    line = queue + line

                                    # Trying to see if both end of line are not connected
                                    pstart = clipped_crown_mesh.vertices[line[0]]
                                    pend = clipped_crown_mesh.vertices[line[-1]]

                                    d = np.linalg.norm(pstart-pend)
                                    # print(d)

                                    open_lines.append(line)
                                    break
                        else:
                            closed_polygons.append(line)
                        break

            except: # FIXME: specifier quelle exception est attendue ici !
                # TODO: retirer les deux lines suivantes
                if self._verbose:
                    print("%u closed polygon\n%u open curve" % (len(closed_polygons), len(open_lines)))

                if self._assert_closed_boundaries:
                    if len(open_lines) > 0:
                        try:
                            from . import mmio
                            mmio.write_VTP('full_debug.vtp', self.source_mesh.vertices, self.source_mesh.faces)
                            mmio.write_VTP('clipped_crown_debug.vtp', clipped_crown_mesh.vertices, clipped_crown_mesh.faces)
                            mmio.write_VTP('crown_debug.vtp', crown_mesh.vertices, crown_mesh.faces)
                        except:
                            pass

                        for line in open_lines:
                            print(line)

                        raise RuntimeError('Open intersection curve found with assert_closed_boundaries option enabled. Files full_debug.vtp, crown_debug.vtp and clipped_crown_debug.vtp written.')

                break

        output = {'clipped_crown_mesh': clipped_crown_mesh,
                  'closed_polygons': closed_polygons,
                  'open_lines': open_lines}

        self.__internals__.update(output)

    def _clip(self):
        """Performs clipping and assemble the clipped mesh."""

        self._clip_crown_by_plane()
        clipped_mesh = self.lower_mesh + self.clipped_crown_mesh
        clipped_mesh.name = '_'.join((self._source_mesh.name, 'clipped'))
        self.__internals__['clipped_mesh'] = clipped_mesh
        return
