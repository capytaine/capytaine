#!/usr/bin/env python
# coding: utf-8
"""
This module concerns mesh data structures.

Based on the code meshmagick from François Rongère (EC Nantes).
"""

from itertools import count

import numpy as np
import vtk

from capytaine.mesh.geometry import Plane
from capytaine.mesh.faces_properties import compute_faces_properties


class Mesh:
    """A class to handle unstructured meshes.

    Parameters
    ----------
    vertices : array_like of shape (nv, 3)
        Array of mesh vertices coordinates. Each line of the array represents one vertex coordinates
    faces : array_like of shape (nf, 4)
        Arrays of mesh connectivities for faces. Each line of the array represents indices of vertices that form the
        face, expressed in counterclockwise order to ensure outward normals description.
    name : str, optional
        The mesh's name. If None, mesh is given an automatic name based on its internal ID.
    """
    _ids = count(0)

    def __init__(self, vertices, faces, name=None):

        self.__internals__ = dict()

        assert np.array(vertices).shape[1] == 3
        assert np.array(faces).shape[1] == 4

        self.vertices = np.array(vertices, dtype=np.float)
        self.faces = np.array(faces, dtype=np.int)
        self._id = next(self._ids)

        if not name:
            self.name = f'mesh_{self._id}'
        else:
            self.name = str(name)

    def __str__(self):
        return self.name

    @property
    def nb_vertices(self):
        """Get the number of vertices in the mesh

        Returns
        -------
        int
        """
        return self.vertices.shape[0]

    @property
    def nb_faces(self):
        """Get the number of faces in the mesh

        Returns
        -------
        int
        """
        return self.faces.shape[0]

    def copy(self, name=None):
        """Get a copy of the current mesh instance.

        Parameters
        ----------
        name : string, optional
            a new name for the new mesh

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

    ##################
    #  Extract face  #
    ##################

    def get_face(self, face_id):
        """Get the face described by its vertices connectivity

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

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """
        Extracts a new mesh from a selection of faces ids

        Parameters
        ----------
        id_faces_to_extract : ndarray
            Indices of faces that have to be extracted
        return_index: bool
            Flag to output old indices

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
        extracted_mesh._verbose = self._verbose

        extracted_mesh.name = 'mesh_extracted_from_%s' % self.name

        if return_index:
            return extracted_mesh, id_v
        else:
            return extracted_mesh

    ######################
    #  Faces properties  #
    ######################

    def _has_faces_properties(self):
        return 'faces_areas' in self.__internals__

    def _remove_faces_properties(self):
        if self._has_faces_properties():
            del self.__internals__['faces_areas']
            del self.__internals__['faces_centers']
            del self.__internals__['faces_normals']
            del self.__internals__['faces_radiuses']
        return

    @property
    def faces_areas(self):
        """Get the array of faces areas of the mesh

        Returns
        -------
        ndarray
        """
        if 'faces_areas' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_areas']

    @property
    def faces_centers(self):
        """Get the array of faces centers of the mesh

        Returns
        -------
        ndarray
        """
        if 'faces_centers' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_centers']

    @property
    def faces_normals(self):
        """Get the array of faces normals of the mesh

        Returns
        -------
        ndarray
        """
        if 'faces_normals' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_normals']

    @property
    def faces_radiuses(self):
        """Get the array of faces radiuses of the mesh.

        Returns
        -------
        ndarray
        """
        if 'faces_radiuses' not in self.__internals__:
            self.__internals__.update(compute_faces_properties(self))
        return self.__internals__['faces_radiuses']

    ###############################
    #  Triangles and quadrangles  #
    ###############################

    def is_triangle(self, face_id):
        """Returns if a face is a triangle

        Parameters
        ----------
        face_id : int
            Face id

        Returns
        -------
        bool
            True if the face with id face_id is a triangle
        """
        assert 0 <= face_id < self.nb_faces
        return self._faces[face_id, 0] == self._faces[face_id, -1]

    def _triangles_quadrangles(self):
        triangle_mask = (self._faces[:, 0] == self._faces[:, -1])
        quadrangles_mask = np.invert(triangle_mask)
        triangles_quadrangles = {'triangles_ids': np.where(triangle_mask)[0],
                                 'quadrangles_ids': np.where(quadrangles_mask)[0]}
        self.__internals__.update(triangles_quadrangles)
        return

    def _has_triangles_quadrangles(self):
        return 'triangles_ids' in self.__internals__

    def _remove_triangles_quadrangles(self):
        if 'triangles_ids' in self.__internals__:
            del self.__internals__['triangles_ids']
            del self.__internals__['quadrangles_ids']
        return

    @property
    def triangles_ids(self):
        """Get the array of ids of triangle shaped faces

        Returns
        -------
        ndarray
        """
        if 'triangles_ids' not in self.__internals__:
            self._triangles_quadrangles()
        return self.__internals__['triangles_ids']

    @property
    def nb_triangles(self):
        """Get the number of triangles in the mesh

        Returns
        -------
        int
        """
        if 'triangles_ids'not in self.__internals__:
            self._triangles_quadrangles()
        return len(self.__internals__['triangles_ids'])

    @property
    def quadrangles_ids(self):
        """Get the array of ids of quadrangle shaped faces

        Returns
        -------
        ndarray
        """
        if 'triangles_ids' not in self.__internals__:
            self._triangles_quadrangles()
        return self.__internals__['quadrangles_ids']

    @property
    def nb_quadrangles(self):
        """Get the number of quadrangles in the mesh

        Returns
        -------
        int
        """
        if 'triangles_ids' not in self.__internals__:
            self._triangles_quadrangles()
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
        This method differs from `axis_aligned_bbox()` by the fact that the bounding box that is returned is squared but have the same center as the AABB
        """
        xmin, xmax, ymin, ymax, zmin, zmax = self.axis_aligned_bbox
        (x0, y0, z0) = np.array([xmin+xmax, ymin+ymax, zmin+zmax]) * 0.5
        d = (np.array([xmax-xmin, ymax-ymin, zmax-zmin]) * 0.5).max()

        return x0-d, x0+d, y0-d, y0+d, z0-d, z0+d

    def show(self, **kwargs):
        self.show_vtk(**kwargs)

    def show_vtk(self):
        """Shows the mesh in the meshmagick viewer"""
        from capytaine.mesh.MMviewer import compute_vtk_polydata, MMViewer

        vtk_polydata = compute_vtk_polydata(self)
        self.viewer = MMViewer()
        self.viewer.add_polydata(vtk_polydata)
        self.viewer.show()
        self.viewer.finalize()

    def show_matplotlib(self, ax=None,
                        normal_vectors=False, scale_normal_vector=None,
                        saveas=None,
                        **kwargs):
        """Poor man's viewer with matplotlib.

        Parameters
        ----------
        ax: matplotlib axis
            The 3d axis in which to plot the mesh. If not provided, create a new one.
        normal_vector: bool
            If True, print normal vector.
        scale_normal_vector: array of shape (nb_faces, )
            Scale separately each of the normal vectors.
        saveas: str
            file path where to save the image

        Other parameters are passed to Poly3DCollection.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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
        if 'facecolors' not in kwargs:
            kwargs['facecolors'] = (0.3, 0.3, 0.3, 0.3)
        if 'edgecolor' not in kwargs:
            kwargs['edgecolor'] = 'k'
        ax.add_collection3d(Poly3DCollection(faces, **kwargs))

        # Plot normal vectors.
        if normal_vectors:
            if scale_normal_vector is not None:
                vectors = (scale_normal_vector * self.faces_normals.T).T
            else:
                vectors = self.faces_normals
            ax.quiver(*zip(*self.faces_centers), *zip(*vectors), length=0.2)

        if default_axis:
            ax.set_xlabel("x")
            ax.set_ylabel("y")

            xmin, xmax, ymin, ymax, zmin, zmax = self.squared_axis_aligned_bbox
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(zmin, zmax)

            if saveas is not None:
                plt.tight_layout()
                plt.savefig(saveas)
            else:
                plt.show()

    ################################
    #  Transformation of the mesh  #
    ################################

    ## Rotate

    def rotate_x(self, thetax):
        """Rotates the mesh around Ox axis.

        Parameters
        ----------
        thetax : float
            Angle (rad)

        Returns
        -------
        ndarray
            The (3x3) rotation matrix that has been applied to rotate the mesh
        """
        return self.rotate([thetax, 0., 0.])

    def rotate_y(self, thetay):
        """Rotates the mesh around Oy axis.

        Parameters
        ----------
        thetay : float
            Angle (rad)

        Returns
        -------
        ndarray
            The (3x3) rotation matrix that has been applied to rotate the mesh
        """
        return self.rotate([0., thetay, 0.])

    def rotate_z(self, thetaz):
        """Rotates the mesh around Oz axis.

        Parameters
        ----------
        thetaz : float
            Angle (rad)

        Returns
        -------
        ndarray
            The (3x3) rotation matrix that has been applied to rotate the mesh
        """
        return self.rotate([0., 0., thetaz])

    def rotate(self, angles):
        """Rotates the mesh in 3D giving the 3 rotation angles that are defined around fixed axes.

        Parameters
        ----------
        angles : array_like
            The 3 angles of the 3D rotation (rad)

        Returns
        -------
        ndarray
            The (3x3) rotation matrix that has been applied to rotate the mesh
        """
        if self.has_surface_integrals():
            self._remove_surface_integrals()
        # TODO: docstring
        # FIXME : code en doublon par rapport a la fonction _rodrigues du module geometry

        angles = np.asarray(angles, dtype=np.float)
        theta = np.linalg.norm(angles)
        if theta == 0.:
            return np.eye(3)

        ctheta = np.cos(theta)
        stheta = np.sin(theta)

        nx, ny, nz = angles/theta
        nxny = nx*ny
        nxnz = nx*nz
        nynz = ny*nz
        nx2 = nx*nx
        ny2 = ny*ny
        nz2 = nz*nz

        rot_matrix = ctheta*np.eye(3) \
            + (1-ctheta) * np.array([[nx2, nxny, nxnz],
                                     [nxny, ny2, nynz],
                                     [nxnz, nynz, nz2]]) \
            + stheta * np.array([[0., -nz, ny],
                                [nz, 0., -nx],
                                [-ny, nx, 0.]])

        # TODO: travailler avec une classe rotation
        self._vertices = np.transpose(np.dot(rot_matrix, self._vertices.copy().T))

        # Updating faces properties if any
        # TODO: use traitlets...
        if self._has_faces_properties():
            # Rotating normals and centers too
            normals = self.__internals__['faces_normals']
            centers = self.__internals__['faces_centers']
            self.__internals__['faces_normals'] = np.transpose(np.dot(rot_matrix, normals.T))
            self.__internals__['faces_centers'] = np.transpose(np.dot(rot_matrix, centers.T))

        # if self.has_surface_integrals():
        #     self._remove_surface_integrals()

        return rot_matrix

    ## Translate

    def translate_x(self, tx):
        """Translates the mesh along the Ox axis.

        Parameters
        ----------
        tx : float
            Distance
        """
        return self.translate([tx, 0.0, 0.0])

    def translate_y(self, ty):
        """Translates the mesh along the Oy axis.

        Parameters
        ----------
        ty : float
            Distance
        """
        return self.translate([0.0, ty, 0.0])

    def translate_z(self, tz):
        """Translates the mesh along the Oz axis.

        Parameters
        ----------
        tz : float
            Distance
        """
        return self.translate([0.0, 0.0, tz])

    def translate(self, t):
        """Translates the mesh in 3D giving the 3 distances along coordinate axes.

        Parameters
        ----------
        t : array_like
            translation vector
        """
        tx, ty, tz = t
        V = self.vertices.copy() # FIXME: why doing a copy ???
        V[:, 0] += tx
        V[:, 1] += ty
        V[:, 2] += tz
        self.vertices = V

        # Updating properties if any
        if self._has_faces_properties():
            centers = self.__internals__['faces_centers']
            centers[:, 0] += tx
            centers[:, 1] += ty
            centers[:, 2] += tz
            self.__internals__['faces_centers'] = centers

        if self.has_surface_integrals():
            self._remove_surface_integrals()

        return

    ## Scale

    def scale(self, alpha):
        """Scales the mesh.

        Parameters
        ----------
        alpha : float
            A positive scaling factor
        """
        assert 0 < alpha

        # TODO: voir pourquoi il est fait une copie ici...
        vertices = self._vertices.copy()
        vertices *= float(alpha)
        self._vertices = vertices

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    def scalex(self, alpha):
        """Scales the mesh along the x axis.

        Parameters
        ----------
        alpha : float
            A positive scaling factor
        """
        assert 0 < alpha

        vertices = self._vertices.copy()
        vertices[:, 0] *= float(alpha)
        self._vertices = vertices

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    def scaley(self, alpha):
        """Scales the mesh along the y axis.

        Parameters
        ----------
        alpha : float
            A positive scaling factor
        """
        assert 0 < alpha

        vertices = self._vertices.copy()
        vertices[:, 1] *= float(alpha)
        self._vertices = vertices

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    def scalez(self, alpha):
        """Scales the mesh along the z axis.

        Parameters
        ----------
        alpha : float
            A positive scaling factor
        """
        assert 0 < alpha

        vertices = self._vertices.copy()
        vertices[:, 2] *= float(alpha)
        self._vertices = vertices

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    ## Other

    def flip_normals(self):
        """Flips every normals of the mesh."""

        faces = self._faces.copy()
        self._faces = np.fliplr(faces)

        if self._has_faces_properties():
            self.__internals__['faces_normals'] *= -1

        if self.has_surface_integrals():
            self._remove_surface_integrals()

        return

    def triangulate_quadrangles(self):
        """Triangulates every quadrangles of the mesh by simple spliting.

        Each quadrangle gives two triangles.

        Note
        ----
        No checking is made on the triangle quality is done.
        """
        # TODO: Ensure the best quality aspect ratio of generated triangles

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

        self.__internals__.clear()

        self._faces = faces

        return faces

    def symmetrize(self, plane):
        """Symmetrize the mesh with respect to a plane.

        Parameters
        ----------
        plane : Plane
            The plane of symmetry
        """
        # Symmetrizing the nodes
        vertices, faces = self._vertices, self._faces

        vertices = np.concatenate((vertices, vertices - 2 * np.outer(np.dot(vertices, plane.normal) - plane.c, plane.normal)))
        faces = np.concatenate((faces, np.fliplr(faces.copy() + self.nb_vertices)))

        self._vertices, self._faces = vertices, faces
        verbose = self.verbose
        self.verbose_off()
        self.merge_duplicates()
        self.verbose = verbose

        self.__internals__.clear()

        return

    def mirror(self, plane):
        """Mirrors the mesh instance with respect to a plane.

        Parameters
        ----------
        plane : Plane
            The mirroring plane
        """
        self._vertices -= 2 * np.outer(np.dot(self._vertices, plane.normal) - plane.c, plane.normal)
        self.flip_normals()
        self.__internals__.clear()
        return

    def get_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Clip the mesh with two horizontal planes."""
        from meshmagick.mesh_clipper import MeshClipper
        if self.vertices[:, 2].min() > free_surface or self.vertices[:, 2].max() < sea_bottom:
            return None  # The mesh has no wet faces.

        elif self.vertices[:, 2].min() > sea_bottom and self.vertices[:, 2].max() < free_surface:
            return self.copy()  # The mesh is completely immersed. Non need for clipping.

        else:
            clipped_mesh = MeshClipper(self,
                                       plane=Plane(normal=(0.0, 0.0, 1.0),
                                                   scalar=free_surface)).clipped_mesh

            if sea_bottom > -np.infty:
                clipped_mesh = MeshClipper(clipped_mesh,
                                           plane=Plane(normal=(0.0, 0.0, -1.0),
                                                       scalar=-sea_bottom)).clipped_mesh

            clipped_mesh.remove_unused_vertices()
            return clipped_mesh

    ####################
    #  Combine meshes  #
    ####################

    def __add__(self, mesh_to_add):
        """Adds two meshes

        Parameters
        ----------
        mesh_to_add : Mesh
            The other mesh instance to add to the current instance

        Returns
        -------
        Mesh
            The composite mesh

        Note
        ----
        This method should not be called as is but it overides the + binary operator for convenience.
        """

        assert isinstance(mesh_to_add, Mesh)
        vertices = np.concatenate((self._vertices, mesh_to_add._vertices), axis=0)
        faces = np.concatenate((self._faces, mesh_to_add._faces + self.nb_vertices), axis=0)
        new_mesh = Mesh(vertices, faces, name='_'.join([self.name, mesh_to_add.name]))
        new_mesh.merge_duplicates()

        return new_mesh

    ##################
    #  Mesh quality  #
    ##################

    def merge_duplicate(self, **kwargs):
        from capytaine.mesh.quality import merge_duplicates
        return merge_duplicates(self, **kwargs)

    def heal_normals(self, **kwargs):
        from capytaine.mesh.quality import heal_normals
        return heal_normals(self, **kwargs)

    def remove_unused_vertices(self, **kwargs):
        from capytaine.mesh.quality import remove_unused_vertices
        return remove_unused_vertices(self, **kwargs)

    def heal_triangles(self, **kwargs):
        from capytaine.mesh.quality import heal_triangles
        return heal_triangles(self, **kwargs)

    def remove_degenerated_faces(self, **kwargs):
        from capytaine.mesh.quality import remove_degenerated_faces
        return remove_degenerated_faces(self, **kwargs)

    def heal_mesh(self):
        """Heals the mesh for different tests available.

        It applies:

        * Unused vertices removal
        * Degenerate faces removal
        * Duplicate vertices merging
        * Triangles healing
        * Normal healing
        """
        if self._has_faces_properties():
            self._remove_faces_properties()
        self.remove_unused_vertices()
        self.remove_degenerated_faces()
        self.merge_duplicates()
        self.heal_triangles()
        self.heal_normals()
        return

