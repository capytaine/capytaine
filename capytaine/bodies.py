#!/usr/bin/env python
# coding: utf-8
"""Floatting bodies to be used in radiation-diffraction problems.

class FloattingBody
class Sphere
class HorizontalCylinder
class OneSidedRectangle
class TwoSidedRectangle
class OpenRectangularParallelepiped
class RectangularParallelepiped
"""

import sys
from itertools import product

import numpy as np

from meshmagick.mesh import Mesh


class FloattingBody(Mesh):
    """A floatting body composed of a mesh (inherited from Meshmagick) and
    several degrees of freedom (dof)."""

    def __init__(self, *args, **kwargs):
        Mesh.__init__(self, *args, **kwargs)
        self.dof = {}

    @staticmethod
    def from_file(filename, file_format):
        from meshmagick.mmio import load_mesh

        vertices, faces = load_mesh(filename, file_format)

        return FloattingBody(vertices, faces, name="mesh_from_"+filename)

    def __add__(self, body_to_add):
        new_body = Mesh.__add__(self, body_to_add)
        new_body.__class__ = FloattingBody
        return new_body

    @property
    def faces_radiuses(self):
        """Get the array of faces radiuses of the mesh

        Returns
        -------
        ndarray
        """
        if 'faces_radiuses' not in self.__internals__:
            self._compute_radiuses()
        return self.__internals__['faces_radiuses']

    def _compute_radiuses(self):
        """Update face radiuses"""
        from numpy.linalg import norm

        faces_radiuses = np.zeros(self.nb_faces, dtype=np.float32)
        for j in range(self.nb_faces): # TODO: optimize by array broadcasting
            faces_radiuses[j] = max(
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 0], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 1], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 2], 0:3]),
                norm(self.faces_centers[j, 0:3] -
                     self.vertices[self.faces[j, 3], 0:3]),
                )

        self.__internals__["faces_radiuses"] = faces_radiuses

    def get_immerged_part(self, depth=np.infty):
        """Use Meshmagick mesh clipper to remove the part of the mesh above the
        free surface and the part of the mesh below the sea bottom.
        """
        from meshmagick.geometry import Plane
        from meshmagick.mesh_clipper import MeshClipper

        clipped_mesh = MeshClipper(self,
                                   plane=Plane(normal=(0.0, 0.0, 1.0),
                                               scalar=0.0)).clipped_mesh

        if depth < np.infty:
            clipped_mesh = MeshClipper(clipped_mesh,
                                       plane=Plane(normal=(0.0, 0.0, -1.0),
                                                   scalar=depth)).clipped_mesh

        clipped_mesh.remove_unused_vertices()

        return FloattingBody(clipped_mesh.vertices, clipped_mesh.faces)

    def show_matplotlib(self, dof=None):
        """Poor man's viewer with matplotlib
        To be deleted when the VTK viewer is fully working with Python 3?"""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        faces = []
        for face in self.faces:
            vertices = []
            for index_vertex in face:
                vertices.append(self.vertices[int(index_vertex), :])
            faces.append(vertices)
        ax.add_collection3d(Poly3DCollection(faces, facecolor=(0.3, 0.3, 0.3, 0.3), edgecolor='k'))

        # Plot normal vectors.
        if dof:
            normals = [self.dof[dof][j] * normal for j, normal in enumerate(self.faces_normals)]
            ax.quiver(*zip(*self.faces_centers), *zip(*normals), length=0.2)

        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(min(self.vertices[:, 0]), max(self.vertices[:, 0]))
        plt.ylim(min(self.vertices[:, 1]), max(self.vertices[:, 1]))
        plt.gca().set_zlim(min(self.vertices[:, 2]), max(self.vertices[:, 2]))
        plt.show()

class Sphere(FloattingBody):
    """Floatting body of the shape of a sphere."""

    def __init__(self, radius=1.0, ntheta=11, nphi=11, z0=0.0, clip_free_surface=False):

        if clip_free_surface:
            if z0 < -radius: # fully immerged
                theta_max = np.pi
            elif z0 < radius:
                theta_max = np.arccos(z0/radius)
            else:
                raise Exception("Sphere out of the water")
        else:
            theta_max = np.pi

        theta = np.linspace(-theta_max, theta_max, ntheta)
        phi = np.linspace(-np.pi/2, np.pi/2, nphi)

        # Nodes
        nodes = np.zeros((ntheta*nphi, 3), dtype=np.float32)

        for i, (t, p) in enumerate(product(theta, phi)):
            # The sign of theta below is a trick to get the correct orientation of the normal vectors...
            x =      radius * np.sin(t) * np.sin(np.sign(t)*p)
            y =      radius * np.sin(t) * np.cos(np.sign(t)*p)
            z = z0 - radius * np.cos(t)
            nodes[i, :] = (x, y, z)

        # Connectivities
        panels = np.zeros(((ntheta-1)*(nphi-1), 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nphi-1))):
            panels[k, :] = (j+i*nphi, j+(i+1)*nphi, j+1+(i+1)*nphi, j+1+i*nphi)

        FloattingBody.__init__(self, nodes, panels)
        # self.merge_duplicates()
        # self.heal_triangles()


class HorizontalCylinder(FloattingBody):
    """Floatting body of the shape of a cylinder of axis Ox."""

    def __init__(self, length=1.0, radius=1.0, nx=11, nr=3, ntheta=11, z0=0.0, clip_free_surface=False):

        if clip_free_surface:
            if z0 < -radius: # fully immerged
                theta_max = np.pi
            elif z0 < radius:
                theta_max = np.arccos(z0/radius)
            else:
                raise Exception("Sphere out of the water")
        else:
            theta_max = np.pi

        theta = np.linspace(-theta_max, theta_max, ntheta)
        X = np.linspace(0.0, length, nx)
        R = np.linspace(0.0, radius, nr)

        # Nodes
        nodes = np.zeros((ntheta*(nx+2*nr), 3), dtype=np.float32)

        for i, (t, x) in enumerate(product(theta, X)):
            y = radius * np.sin(t)
            z = z0 - radius * np.cos(t)
            nodes[i, :] = (x, y, z)

        for i, (x, r, t) in enumerate(product([0, length], R, theta)):
            y = r * np.sin(t)
            z = z0 - r * np.cos(t)
            nodes[ntheta*nx+i, :] = (x, y, z)

        # Connectivities
        npanels = (ntheta-1)*((nx-1)+2*(nr-1))
        panels = np.zeros((npanels, 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nx-1))):
            panels[k, :] = (j+i*nx, j+(i+1)*nx, j+1+(i+1)*nx, j+1+i*nx)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*nx, ntheta*nx+ntheta-1))):
            panels[(ntheta-1)*(nx-1)+k, :] = (j+i*ntheta, j+1+i*ntheta, j+1+(i+1)*ntheta, j+(i+1)*ntheta)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            panels[(ntheta-1)*((nx-1)+(nr-1))+k, :] = (j+i*ntheta, j+(i+1)*ntheta, j+1+(i+1)*ntheta, j+1+i*ntheta)

        FloattingBody.__init__(self, nodes, panels)
        self.merge_duplicates()
        self.heal_triangles()


class OneSidedRectangle(FloattingBody):
    """Rectangular panel with cartesian mesh."""

    def __init__(self, height=2.0, length=10.0, nh=5, nl=5, z0=0.0):

        X = np.linspace(-length/2, length/2, nl)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nl*nh, 3), dtype=np.float32)
        panels = np.zeros(((nl-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nl-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)

        FloattingBody.__init__(self, nodes, panels)


class TwoSidedRectangle(FloattingBody):
    """Rectangular panel with cartesian mesh.
    Each face is defined twice with two opposite normal vectors."""

    def __init__(self, height=2.0, length=10.0, nh=5, nl=5, z0=0.0):

        X = np.linspace(-length/2, length/2, nl)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nl*nh, 3), dtype=np.float32)
        panels = np.zeros((2*(nl-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nl-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)
            panels[(nl-1)*(nh-1)+k, :] = (j+i*nh, j+(i+1)*nh, j+1+(i+1)*nh, j+1+i*nh)

        FloattingBody.__init__(self, nodes, panels)


class OpenRectangularParallelepiped(FloattingBody):
    """Four panels forming a parallelepiped without top nor bottom."""

    def __init__(self, height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=3, z0=0.0):
        front = OneSidedRectangle(height=height, length=width, nh=nh, nl=nw, z0=z0)
        back = front.copy()

        front.translate_y(thickness/2)
        back.rotate_z(np.pi)
        back.translate_y(-thickness/2)

        side = OneSidedRectangle(height=height, length=thickness, nh=nh, nl=nth, z0=z0)
        other_side = side.copy()

        side.rotate_z(np.pi/2)
        side.translate_x(-width/2)
        other_side.rotate_z(-np.pi/2)
        other_side.translate_x(width/2)

        combine = front + side + other_side + back
        combine.merge_duplicates()
        combine.heal_triangles()

        FloattingBody.__init__(self, combine.vertices, combine.faces)


class RectangularParallelepiped(FloattingBody):
    """Six panels forming a complete parallelepiped."""

    def __init__(self, height=10.0, width=10.0, thickness=2.0, nh=5, nl=5, nth=3):
        sides = OpenRectangularParallelepiped(height=height, width=width, thickness=thickness, nh=nh, nl=nl, nth=nth)
        top = OneSidedRectangle(height=thickness, length=width, nh=nth, nl=nl)
        bottom = top.copy()

        top.rotate_x(np.pi/2)
        top.translate_y(thickness/2)
        top.translate_z(height)
        bottom.rotate_x(-np.pi/2)
        bottom.translate_y(-thickness/2)

        combine = sides + top + bottom
        combine.merge_duplicates()
        combine.heal_triangles()

        FloattingBody.__init__(self, combine.vertices, combine.faces)
