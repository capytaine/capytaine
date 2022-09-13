#!/usr/bin/env python
# coding: utf-8
""" This module contains a class to describe the 2D mesh of the surface of a body in a 3D space.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from itertools import count

import numpy as np
from numpy.linalg import norm

from capytaine.meshes.geometry import Abstract3DObject, Plane, inplace_transformation
from capytaine.meshes.properties import compute_faces_properties, compute_connectivity
from capytaine.meshes.surface_integrals import compute_faces_integrals
from capytaine.meshes.quality import (merge_duplicates, heal_normals, remove_unused_vertices,
                                      heal_triangles, remove_degenerated_faces)
from capytaine.tools.optional_imports import import_optional_dependency

LOG = logging.getLogger(__name__)


class Mesh(Abstract3DObject):
    """A class to handle unstructured 2D meshes in a 3D space.

    Parameters
    ----------
    vertices : array_like of shape (nv, 3)
        Array of mesh vertices coordinates.Each line of the array represents one vertex
        coordinates
    faces : array_like of shape (nf, 4)
        Arrays of mesh connectivities for faces. Each line of the array represents indices of
        vertices that form the face, expressed in counterclockwise order to ensure outward normals
        description.
    name : str, optional
        The name of the mesh. If None, the mesh is given an automatic name based on its internal ID.
    """

    _ids = count(0)  # A counter for automatic naming of new meshes.

    def __init__(self, vertices=None, faces=None, name=None):

        if vertices is None or len(vertices) == 0:
            vertices = np.zeros((0, 3))

        if faces is None or len(faces) == 0:
            faces = np.zeros((0, 4))

        if name is None:
            self.name = f'mesh_{next(Mesh._ids)}'
        else:
            self.name = str(name)

        self.__internals__ = dict()
        self.vertices = vertices  # Not a direct assignment, goes through the setter method below.
        self.faces = faces  # Not a direct assignment, goes through the setter method below.

        LOG.debug(f"New mesh: {repr(self)}")

    def __str__(self):
        return self.name

    def __repr__(self):
        return (f"{self.__class__.__name__}(nb_vertices={self.nb_vertices}, "
                f"nb_faces={self.nb_faces}, name={self.name})")

    @property
    def nb_vertices(self) -> int:
        """Get the number of vertices in the mesh."""
        return self._vertices.shape[0]

    @property
    def vertices(self) -> np.ndarray:
        """Get the vertices array coordinate of the mesh."""
        return self._vertices

    @vertices.setter
    def vertices(self, value) -> None:
        self._vertices = np.array(value, dtype=float)
        assert self._vertices.shape[1] == 3, \
            "Vertices of a mesh should be provided as a sequence of 3-ple."
        self.__internals__.clear()

    @property
    def nb_faces(self) -> int:
        """Get the number of faces in the mesh."""
        return self._faces.shape[0]

    @property
    def faces(self) -> np.ndarray:
        """Get the faces connectivity array of the mesh."""
        return self._faces

    @faces.setter
    def faces(self, faces):
        faces = np.array(faces, dtype=int)
        assert np.all(faces >= 0), \
            "Faces of a mesh should be provided as positive integers (ids of vertices)"
        assert faces.shape[1] == 4, \
            "Faces of a mesh should be provided as a sequence of 4-ple."
        assert len(faces) == 0 or faces.max()+1 <= self.nb_vertices, \
            "The array of faces should only reference vertices that are in the mesh."
        self._faces = faces
        self.__internals__.clear()

    def copy(self, name=None) -> 'Mesh':
        """Get a copy of the current mesh instance.

        Parameters
        ----------
        name : string, optional
            a name for the new mesh

        Returns
        -------
        Mesh
            mesh instance copy
        """
        from copy import deepcopy
        new_mesh = deepcopy(self)
        if name is not None:
            new_mesh.name = name
        return new_mesh

    def merged(self):
        """Dummy method to be generalized for collections of meshes."""
        return self

    def tree_view(self, **kwargs):
        """Dummy method to be generalized for collections of meshes."""
        return self.name

    def to_meshmagick(self):
        """Convert the Mesh object as a Mesh object from meshmagick.
        Mostly for debugging."""
        from meshmagick.mesh import Mesh
        meshmagick_mesh = Mesh(self.vertices, self.faces, name=self.name)
        meshmagick_mesh.heal_mesh()
        return meshmagick_mesh

    ##################
    #  Extract face  #
    ##################

    def get_face(self, face_id):
        """Get the face described by its vertices connectivity.

        Parameters
        ----------
        face_id : int
            Face id

        Returns
        -------
        ndarray
            If the face is a triangle, the array has 3 components, else it has 4 (quadrangle)
        """
        if self.is_triangle(face_id):
            return self._faces[face_id, :3]
        else:
            return self._faces[face_id]

    def extract_one_face(self, id_face):
        vertices = self.vertices[self.faces[id_face, :], :]
        mesh = SingleFace(vertices)

        for prop in self.__internals__:
            if prop[:4] == "face":
                mesh.__internals__[prop] = self.__internals__[prop][[id_face]]

        return mesh

    def extract_faces(self, id_faces_to_extract, return_index=False, name=None):
        """
        Extracts a new mesh from a selection of faces ids

        Parameters
        ----------
        id_faces_to_extract : ndarray
            Indices of faces that have to be extracted
        return_index: bool, optional
            Flag to output old indices
        name: string, optional
            Name for the new mesh

        Returns
        -------
        Mesh
            A new Mesh instance composed of the extracted faces
        """
        nv = self.nb_vertices

        # Determination of the vertices to keep
        vertices_mask = np.zeros(nv, dtype=bool)
        vertices_mask[self._faces[id_faces_to_extract].flatten()] = True
        id_v = np.arange(nv)[vertices_mask]

        # Building up the vertex array
        v_extracted = self._vertices[id_v]
        new_id__v = np.arange(nv)
        new_id__v[id_v] = np.arange(len(id_v))

        faces_extracted = self._faces[id_faces_to_extract]
        faces_extracted = new_id__v[faces_extracted.flatten()].reshape((len(id_faces_to_extract), 4))

        extracted_mesh = Mesh(v_extracted, faces_extracted)

        for prop in self.__internals__:
            if prop[:4] == "face":
                extracted_mesh.__internals__[prop] = self.__internals__[prop][id_faces_to_extract]

        if name is None:
            if self.name is not None and self.name.startswith("mesh_extracted_from_"):
                extracted_mesh.name = self.name
            else:
                extracted_mesh.name = f"mesh_extracted_from_{self.name}"
        else:
            extracted_mesh.name = name

        if return_index:
            return extracted_mesh, id_v
        else:
            return extracted_mesh

    def sliced_by_plane(self, plane: Plane):
        from capytaine.meshes.collections import CollectionOfMeshes
        faces_ids_on_one_side = np.where(plane.distance_to_point(self.faces_centers) < 0)[0]
        if len(faces_ids_on_one_side) == 0 or len(faces_ids_on_one_side) == self.nb_faces:
            return self.copy()
        else:
            mesh_part_1 = self.extract_faces(faces_ids_on_one_side)
            mesh_part_2 = self.extract_faces(list(set(range(self.nb_faces)) - set(faces_ids_on_one_side)))
            return CollectionOfMeshes([mesh_part_1, mesh_part_2],
                                      name=f"{self.name}_splitted_by_{plane}")


    #####################
    #  Mean and radius  #
    #####################

    @property
    def center_of_mass_of_nodes(self):
        """(Non-weighted) center of mass of the nodes of the mesh."""
        if 'center_of_mass_of_nodes' not in self.__internals__:
            center_of_mass_of_nodes = np.mean(self.vertices, axis=0)
            self.__internals__['center_of_mass_of_nodes'] = center_of_mass_of_nodes
            return center_of_mass_of_nodes
        return self.__internals__['center_of_mass_of_nodes']

    @property
    def diameter_of_nodes(self):
        """Maximum distance between two nodes of the mesh."""
        if 'diameter_of_nodes' not in self.__internals__:
            diameter_of_nodes = 2*np.max(
                np.linalg.norm(self.vertices - self.center_of_mass_of_nodes, axis=-1)
            )
            self.__internals__['diameter_of_nodes'] = diameter_of_nodes
            return diameter_of_nodes
        return self.__internals__['diameter_of_nodes']

    ######################
    #  Faces properties  #
    ######################

    @property
    def faces_areas(self) -> np.ndarray:
        """Get the array of faces areas of the mesh."""
        if 'faces_areas' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_areas']

    @property
    def faces_centers(self) -> np.ndarray:
        """Get the array of faces centers of the mesh."""
        if 'faces_centers' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_centers']

    @property
    def faces_normals(self) -> np.ndarray:
        """Get the array of faces normals of the mesh."""
        if 'faces_normals' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_normals']

    @property
    def faces_radiuses(self) -> np.ndarray:
        """Get the array of faces radiuses of the mesh."""
        if 'faces_radiuses' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_radiuses']

    @property
    def quadrature_points(self):
        if 'quadrature' in self.__internals__:
            return self.__internals__['quadrature']
        else:
            # Default: first order quadrature
            return (
                self.faces_centers.reshape((self.nb_faces, 1, 3)),  # Points
                self.faces_areas.reshape((self.nb_faces, 1))        # Weights
            )

    @property
    def quadrature_method(self):
        if 'quadrature_method' in self.__internals__:
            return self.__internals__['quadrature_method']
        else:
            return None

    def compute_quadrature(self, method):
        quadpy = import_optional_dependency("quadpy")
        transform = quadpy.c2.transform
        get_detJ = quadpy.cn._helpers.get_detJ

        if method is None:
            # No quadrature (i.e. default first order quadrature)
            if 'quadrature' in self.__internals__:
                del self.__internals__['quadrature']
                del self.__internals__['quadrature_method']
            else:
                pass

        elif isinstance(method, quadpy.c2._helpers.C2Scheme):
            assert method.points.shape[0] == method.dim == 2
            nb_points = method.points.shape[1]
            points = np.empty((self.nb_faces, nb_points, 3))
            weights = np.empty((self.nb_faces, nb_points))

            self.heal_triangles()

            for i_face in range(self.nb_faces):
                # Define a local frame (Oxyz) such that
                # * the corner A of the quadrilateral panel is the origin of the local frame
                # * the edge AB of the quadrilateral panel is along the local x-axis,
                # * the quadrilateral panel is within the local xy-plane (that is, its normal is along the local z-axis).
                # Hence, the corners of the panels all have 0 as z-coordinate in the local frame.

                # Coordinates in global frame
                global_A, global_B, global_C, global_D = self.vertices[self.faces[i_face, :], :]
                n = self.faces_normals[i_face, :]

                ex = (global_B-global_A)/norm(global_B-global_A)  # unit vector of the local x-axis
                ez = n/norm(n)                                    # unit vector of the local z-axis
                ey = np.cross(ex, ez)                             # unit vector of the local y-axis, such that the basis is orthonormal

                R = np.array([ex, ey, ez])
                local_A = np.zeros((3,))             # coordinates of A in local frame, should be zero by construction
                local_B = R @ (global_B - global_A)  # coordinates of B in local frame
                local_C = R @ (global_C - global_A)  # coordinates of C in local frame
                local_D = R @ (global_D - global_A)  # coordinates of D in local frame

                local_quadrilateral = np.array([[local_A, local_D], [local_B, local_C]])[:, :, :-1]
                # Removing last index in last dimension because not interested in z-coordinate which is 0.

                local_quadpoints = transform(method.points, local_quadrilateral)

                local_quadpoints_in_3d = np.concatenate([local_quadpoints, np.zeros((nb_points, 1))], axis=1)
                global_quadpoints = np.array([R.T @ p for p in local_quadpoints_in_3d]) + global_A
                points[i_face, :, :] = global_quadpoints

                weights[i_face, :] = method.weights * 4 * np.abs(get_detJ(method.points, local_quadrilateral))

            self.__internals__['quadrature'] = (points, weights)
            self.__internals__['quadrature_method'] = method

        else:
            raise NotImplementedError

    ###############################
    #  Triangles and quadrangles  #
    ###############################

    def is_triangle(self, face_id) -> bool:
        """Returns if a face is a triangle

        Parameters
        ----------
        face_id : int
            Face id
        """
        assert 0 <= face_id < self.nb_faces
        return self._faces[face_id, 0] == self._faces[face_id, -1]

    def _compute_triangles_quadrangles(self):
        triangle_mask = (self._faces[:, 0] == self._faces[:, -1])
        quadrangles_mask = np.invert(triangle_mask)
        triangles_quadrangles = {'triangles_ids': np.where(triangle_mask)[0],
                                 'quadrangles_ids': np.where(quadrangles_mask)[0]}
        self.__internals__.update(triangles_quadrangles)

    @property
    def triangles_ids(self) -> np.ndarray:
        """Get the array of ids of triangle shaped faces."""
        if 'triangles_ids' not in self.__internals__:
            self._compute_triangles_quadrangles()
        return self.__internals__['triangles_ids']

    @property
    def nb_triangles(self) -> int:
        """Get the number of triangles in the mesh."""
        if 'triangles_ids'not in self.__internals__:
            self._compute_triangles_quadrangles()
        return len(self.__internals__['triangles_ids'])

    @property
    def quadrangles_ids(self) -> np.ndarray:
        """Get the array of ids of quadrangle shaped faces."""
        if 'triangles_ids' not in self.__internals__:
            self._compute_triangles_quadrangles()
        return self.__internals__['quadrangles_ids']

    @property
    def nb_quadrangles(self) -> int:
        """Get the number of quadrangles in the mesh."""
        if 'triangles_ids' not in self.__internals__:
            self._compute_triangles_quadrangles()
        return len(self.__internals__['quadrangles_ids'])

    #############
    #  Display  #
    #############

    @property
    def axis_aligned_bbox(self):
        """Get the axis aligned bounding box of the mesh.

        Returns
        -------
        tuple
            (xmin, xmax, ymin, ymax, zmin, zmax)
        """
        if self.nb_vertices > 0:
            x, y, z = self._vertices.T
            return (x.min(), x.max(),
                    y.min(), y.max(),
                    z.min(), z.max())
        else:
            return tuple(np.zeros(6))

    @property
    def squared_axis_aligned_bbox(self):
        """Get a squared axis aligned bounding box of the mesh.

        Returns
        -------
        tuple
            (xmin, xmax, ymin, ymax, zmin, zmax)

        Note
        ----
        This method differs from `axis_aligned_bbox()` by the fact that
        the bounding box that is returned is squared but have the same center as the `axis_aligned_bbox()`.
        """
        xmin, xmax, ymin, ymax, zmin, zmax = self.axis_aligned_bbox
        (x0, y0, z0) = np.array([xmin+xmax, ymin+ymax, zmin+zmax]) * 0.5
        d = (np.array([xmax-xmin, ymax-ymin, zmax-zmin]) * 0.5).max()

        return x0-d, x0+d, y0-d, y0+d, z0-d, z0+d

    def show(self, **kwargs):
        self.show_vtk(**kwargs)

    def show_vtk(self, **kwargs):
        """Shows the mesh in the vtk viewer"""
        from capytaine.ui.vtk.mesh_viewer import MeshViewer

        viewer = MeshViewer()
        viewer.add_mesh(self, **kwargs)
        viewer.show()
        viewer.finalize()

    def show_matplotlib(self, ax=None,
                        normal_vectors=False, scale_normal_vector=None,
                        saveas=None, color_field=None, cmap=None,
                        cbar_label=None,
                        **kwargs):
        """Poor man's viewer with matplotlib.

        Parameters
        ----------
        ax: matplotlib axis
            The 3d axis in which to plot the mesh. If not provided, create a new one.
        normal_vectors: bool
            If True, print normal vector.
        scale_normal_vector: array of shape (nb_faces, )
            Scale separately each of the normal vectors.
        saveas: str
            File path where to save the image.
        color_field: array of shape (nb_faces, )
            Scalar field to be plot on the mesh (optional).
        cmap: matplotlib colormap
            Colormap to use for field plotting.
        cbar_label: string
            Label for colormap

        Other parameters are passed to Poly3DCollection.
        """
        matplotlib = import_optional_dependency("matplotlib")
        plt = matplotlib.pyplot
        cm = matplotlib.cm

        mpl_toolkits = import_optional_dependency("mpl_toolkits", package_name="matplotlib")
        Poly3DCollection = mpl_toolkits.mplot3d.art3d.Poly3DCollection

        default_axis = ax is None
        if default_axis:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        faces = []
        for face in self.faces:
            vertices = []
            for index_vertex in face:
                vertices.append(self.vertices[int(index_vertex), :])
            faces.append(vertices)

        if color_field is None:
            if 'facecolors' not in kwargs:
                kwargs['facecolors'] = "yellow"
        else:
            if cmap is None:
                cmap = cm.get_cmap('coolwarm')
            m = cm.ScalarMappable(cmap=cmap)
            m.set_array([min(color_field), max(color_field)])
            m.set_clim(vmin=min(color_field), vmax=max(color_field))
            colors = m.to_rgba(color_field)
            kwargs['facecolors'] = colors
        if 'edgecolor' not in kwargs:
            kwargs['edgecolor'] = 'k'

        ax.add_collection3d(Poly3DCollection(faces, **kwargs))

        if color_field is not None:
            cbar = plt.colorbar(m)
            if cbar_label is not None:
                cbar.set_label(cbar_label)



        # Plot normal vectors.
        if normal_vectors:
            if scale_normal_vector is not None:
                vectors = (scale_normal_vector * self.faces_normals.T).T
            else:
                vectors = self.faces_normals
            ax.quiver(*zip(*self.faces_centers), *zip(*vectors), length=0.2)


        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        xmin, xmax, ymin, ymax, zmin, zmax = self.squared_axis_aligned_bbox
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(zmin, zmax)

        if default_axis:
            if saveas is not None:
                plt.tight_layout()
                plt.savefig(saveas)
            else:
                plt.show()

    ################################
    #  Transformation of the mesh  #
    ################################

    @inplace_transformation
    def translate(self, vector) -> 'Mesh':
        """Translates the mesh in 3D giving the 3 distances along coordinate axes.

        Parameters
        ----------
        vector : array_like
            translation vector
        """
        vector = np.asarray(vector, dtype=float)
        assert vector.shape == (3,), "The translation vector should be given as a 3-ple of values."

        self.vertices += vector

        return self

    @inplace_transformation
    def rotate(self, axis, angle) -> 'Mesh':
        """Rotate the mesh of a given angle around an axis.

        Parameters
        ----------
        axis : Axis
        angle : float
        """

        self._vertices = axis.rotate_points(self._vertices, angle)

        return self

    # OTHER
    @inplace_transformation
    def flip_normals(self) -> 'Mesh':
        """Flips every normals of the mesh."""

        self._faces = np.fliplr(self._faces)

        return self

    @inplace_transformation
    def mirror(self, plane) -> 'Mesh':
        """Flip the mesh with respect to a plane.

        Parameters
        ----------
        plane : Plane
            The mirroring plane
        """
        self.vertices -= 2 * np.outer(np.dot(self.vertices, plane.normal) - plane.c, plane.normal)
        self.flip_normals()
        return self

    @inplace_transformation
    def clip(self, plane) -> 'Mesh':
        from capytaine.meshes.clipper import clip
        clipped_self = clip(self, plane=plane)
        self.vertices = clipped_self.vertices
        self.faces = clipped_self.faces
        self._clipping_data = clipped_self._clipping_data
        return self

    def clipped(self, plane, **kwargs) -> 'Mesh':
        # Same API as for the other transformations
        return self.clip(plane, inplace=False, **kwargs)

    def symmetrized(self, plane):
        from capytaine.meshes.symmetric import ReflectionSymmetricMesh
        half = self.clipped(plane, name=f"{self.name}_half")
        return ReflectionSymmetricMesh(half, plane=plane, name=f"symmetrized_of_{self.name}")

    @inplace_transformation
    def keep_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Clip the mesh with two horizontal planes corresponding
        with the free surface and the sea bottom."""
        self.clip(Plane(normal=(0, 0, 1), point=(0, 0, free_surface)))
        if sea_bottom > -np.infty:
            self.clip(Plane(normal=(0, 0, -1), point=(0, 0, sea_bottom)))
        return self

    @inplace_transformation
    def triangulate_quadrangles(self) -> 'Mesh':
        """Triangulates every quadrangles of the mesh by simple splitting.
        Each quadrangle gives two triangles.

        Note
        ----
        No checking on the triangle quality is done.
        """
        # Defining both triangles id lists to be generated from quadrangles
        t1 = (0, 1, 2)
        t2 = (0, 2, 3)

        faces = self._faces

        # Triangulation
        new_faces = faces[self.quadrangles_ids].copy()
        new_faces[:, :3] = new_faces[:, t1]
        new_faces[:, -1] = new_faces[:, 0]

        faces[self.quadrangles_ids, :3] = faces[:, t2][self.quadrangles_ids]
        faces[self.quadrangles_ids, -1] = faces[self.quadrangles_ids, 0]

        faces = np.concatenate((faces, new_faces))

        LOG.info('\nTriangulating quadrangles')
        if self.nb_quadrangles != 0:
            LOG.info('\t-->{:d} quadrangles have been split in triangles'.format(self.nb_quadrangles))

        self._faces = faces

        return self

    ####################
    #  Combine meshes  #
    ####################

    def join_meshes(*meshes, name=None):
        from capytaine.meshes.collections import CollectionOfMeshes
        return CollectionOfMeshes(meshes, name=name).merged()

    def __add__(self, mesh_to_add) -> 'Mesh':
        return self.join_meshes(mesh_to_add)

    ####################
    #  Compare meshes  #
    ####################
    # The objective is to write a mesh as a set of faces in order to check for equality or to
    # compute differences of meshes. Each face can be represented as a 4x3 array (4 triplets of
    # coordinates).
    # However, it is tricky on several aspects:
    #   * The builtin set class compares the hashes of its objects. Since numpy ndarray are not
    #   hashable, the 4x3 array can be transformed into a tuple of tuples (which is hashable).
    #   * Two faces with different numbering of the faces (but the same ordering) are recorded as
    #   different.
    #   * Two vertices equal up to machine precision can be recorded as different, due to rounding
    #   errors.
    #
    # A possible solution is to define a Face class and a Vertex class with the appropriate __eq__
    # and __hash__.
    #
    # The current implementation below is a rough draft.
    # However, the equality shall still be use for testing.

    def as_set_of_faces(self):
        return frozenset(frozenset(tuple(vertex) for vertex in face) for face in self.vertices[self.faces])

    @staticmethod
    def from_set_of_faces(set_of_faces):
        faces = []
        vertices = []
        for face in set_of_faces:
            ids_of_vertices_in_face = []

            for vertex in face:
                if vertex not in vertices:
                    i = len(vertices)
                    vertices.append(vertex)
                else:
                    i = vertices.index(vertex)
                ids_of_vertices_in_face.append(i)

            if len(ids_of_vertices_in_face) == 3:
                # Add a fourth node identical to the first one
                ids_of_vertices_in_face.append(ids_of_vertices_in_face[0])

            faces.append(ids_of_vertices_in_face)
        return Mesh(vertices=vertices, faces=faces)

    def __eq__(self, other):
        if not isinstance(other, Mesh):
            return NotImplemented
        else:
            return self.as_set_of_faces() == other.as_set_of_faces()

    def __hash__(self):
        if 'hash' not in self.__internals__:
            self.__internals__['hash'] = hash(self.as_set_of_faces())
        return self.__internals__['hash']

    ##################
    #  Mesh quality  #
    ##################

    def merge_duplicates(self, **kwargs):
        return merge_duplicates(self, **kwargs)

    def heal_normals(self, **kwargs):
        return heal_normals(self, **kwargs)

    def remove_unused_vertices(self, **kwargs):
        return remove_unused_vertices(self, **kwargs)

    def heal_triangles(self, **kwargs):
        return heal_triangles(self, **kwargs)

    def remove_degenerated_faces(self, **kwargs):
        return remove_degenerated_faces(self, **kwargs)

    @inplace_transformation
    def heal_mesh(self, closed_mesh=True):
        """Heals the mesh for different tests available.

        It applies:

        * Unused vertices removal
        * Degenerate faces removal
        * Duplicate vertices merging
        * Triangles healing
        * Normal healing
        """
        self.remove_unused_vertices()
        self.remove_degenerated_faces()
        self.merge_duplicates()
        self.heal_triangles()
        if closed_mesh:
            self.heal_normals()
        return self

    #################
    #  Edges stats  #
    #################

    def _edges_stats(self):
        """Computes the min, max, and mean of the mesh's edge length"""
        vertices = self.vertices[self.faces]
        edge_length = np.zeros((self.nb_faces, 4), dtype=float)
        for i in range(4):
            edge = vertices[:, i, :] - vertices[:, i-1, :]
            edge_length[:, i] = np.sqrt(np.einsum('ij, ij -> i', edge, edge))

        return edge_length.min(), edge_length.max(), edge_length.mean()

    @property
    def min_edge_length(self) -> float:
        """The mesh's minimum edge length"""
        return self._edges_stats()[0]

    @property
    def max_edge_length(self) -> float:
        """The mesh's maximum edge length"""
        return self._edges_stats()[1]

    @property
    def mean_edge_length(self) -> float:
        """The mesh's mean edge length"""
        return self._edges_stats()[2]

    #######################
    #  Surface integrals  #
    #######################

    def get_surface_integrals(self) -> np.ndarray:
        """Get the mesh surface integrals."""
        if 'surface_integrals' not in self.__internals__:
            self.__internals__['surface_integrals'] = compute_faces_integrals(self)
        return self.__internals__['surface_integrals']

    @property
    def volume(self) -> float:
        """Get the mesh enclosed volume."""
        normals = self.faces_normals
        sigma_0_2 = self.get_surface_integrals()[:3]

        return (normals.T * sigma_0_2).sum() / 3.

    ####################
    #  Connectivities  #
    ####################

    @property
    def vv(self) -> dict:
        """Get the vertex / vertex connectivity dictionary."""
        if 'v_v' not in self.__internals__:
            self.__internals__.update(compute_connectivity(self))
        return self.__internals__['v_v']

    @property
    def vf(self) -> dict:
        """Get the vertex / faces connectivity dictionary."""
        if 'v_f' not in self.__internals__:
            self.__internals__.update(compute_connectivity(self))
        return self.__internals__['v_f']

    @property
    def ff(self) -> dict:
        """Get the face / faces connectivity dictionary."""
        if 'f_f' not in self.__internals__:
            self.__internals__.update(compute_connectivity(self))
        return self.__internals__['f_f']

    @property
    def boundaries(self) -> list:
        """Get a list that stores lists of boundary connected vertices."""
        if 'boundaries' not in self.__internals__:
            self.__internals__.update(compute_connectivity(self))
        return self.__internals__['boundaries']

    @property
    def nb_boundaries(self) -> int:
        """Get the number of boundaries in the mesh."""
        if 'boundaries' not in self.__internals__:
            self.__internals__.update(compute_connectivity(self))
        return len(self.__internals__['boundaries'])


class SingleFace(Mesh):
    """A view on a single face of a mesh.
    To be used for ACA."""

    _faces = np.arange(4).reshape((1, 4))
    name = "some single face"

    def __init__(self, vertices=None):
        self._vertices = vertices
        self.__internals__ = dict()
