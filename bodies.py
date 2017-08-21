#!/usr/bin/env python
# coding: utf-8

import sys
from itertools import product

import numpy as np
from numpy.linalg import norm

sys.path.append("/home/ancellin/meshmagick/")
from meshmagick.mesh import Mesh, Plane
from meshmagick.mesh_clipper import MeshClipper


class FloattingBody(Mesh):
    def __init__(self, *args, **kwargs):
        Mesh.__init__(self, *args, **kwargs)
        self.dof = {}

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
        faces_radiuses = np.zeros(self.nb_faces, dtype=np.float32)
        for j in range(self.nb_faces): # TODO: optimize by array broadcasting
            faces_radiuses[j] = max(
                    norm(self.faces_centers[j, 0:3] - self.vertices[self.faces[j, 0], 0:3]),
                    norm(self.faces_centers[j, 0:3] - self.vertices[self.faces[j, 1], 0:3]),
                    norm(self.faces_centers[j, 0:3] - self.vertices[self.faces[j, 2], 0:3]),
                    norm(self.faces_centers[j, 0:3] - self.vertices[self.faces[j, 3], 0:3]),
                    )

        self.__internals__["faces_radiuses"] = faces_radiuses


    def keep_only_immerged_part(self, depth=np.infty):
        """Use Meshmagick mesh clipper to remove the part of the mesh above the
        free surface and he part of the mesh below the sea bottom.
        """
        clipped_mesh = MeshClipper(self, plane=Plane(normal=(0.0, 0.0, 1.0), scalar=0.0)).clipped_mesh

        if depth < np.infty:
            clipped_mesh = MeshClipper(clipped_mesh, plane=Plane(normal=(0.0, 0.0, -1.0), scalar=depth)).clipped_mesh

        self.vertices = clipped_mesh.vertices
        self.faces    = clipped_mesh.faces


class Sphere(FloattingBody):
    def __init__(self, radius=1.0, ntheta=11, nphi=11, z0=0.0):

        theta = np.linspace(-np.pi, np.pi, ntheta)
        phi = np.linspace(-np.pi/2, np.pi/2, nphi)

        # Nodes
        nodes  = np.zeros((ntheta*nphi, 3), dtype=np.float32)

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
        self.merge_duplicates()
        self.heal_triangles()


class HorizontalCylinder(FloattingBody):
    def __init__(self, length=1.0, radius=1.0, nx=11, nr=3, ntheta=11, z0=0.0):

        X = np.linspace(0.0, length, nx)
        R = np.linspace(0.0, radius, nr)
        theta = np.linspace(-np.pi, np.pi, ntheta)

        # Nodes
        nnodes = ntheta*(nx+2*nr)
        nodes  = np.zeros((nnodes, 3), dtype=np.float32)

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
            panels[(ntheta-1)*(nx-1)+k, :] = (j+i*ntheta,  j+1+i*ntheta, j+1+(i+1)*ntheta, j+(i+1)*ntheta)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            panels[(ntheta-1)*((nx-1)+(nr-1))+k, :] = (j+i*ntheta, j+(i+1)*ntheta, j+1+(i+1)*ntheta, j+1+i*ntheta)

        FloattingBody.__init__(self, nodes, panels)
        self.merge_duplicates()
        self.heal_triangles()


class OneSidedRectangle(FloattingBody):
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
    def __init__(self, height=2.0, length=10.0, nh=5, nl=5, z0=0.0):

        X = np.linspace(-length/2, length/2, nl)
        Z = np.linspace(z0, z0+height, nh)

        nodes   = np.zeros((nl*nh, 3), dtype=np.float32)
        panels  = np.zeros((2*(nl-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nl-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)
            panels[(nl-1)*(nh-1)+k, :] = (j+i*nh, j+(i+1)*nh, j+1+(i+1)*nh,  j+1+i*nh, )

        FloattingBody.__init__(self, nodes, panels)


class OpenRectangularParallelepiped(FloattingBody):
    def __init__(self, height=10.0, length=10.0, thickness=2.0, nh=5, nl=5, nth=3, z0=0.0):
        front = OneSidedRectangle(height=height, length=length, nh=nh, nl=nl, z0=z0)
        back = front.copy()

        front.translate_y(thickness/2)
        back.rotate_z(np.pi)
        back.translate_y(-thickness/2)

        side  = OneSidedRectangle(height=height, length=thickness, nh=nh, nl=nth, z0=z0)
        other_side = side.copy()

        side.rotate_z(np.pi/2)
        side.translate_x(-length/2)
        other_side.rotate_z(-np.pi/2)
        other_side.translate_x(length/2)

        combine = front + side + other_side + back
        combine.merge_duplicates()
        combine.heal_triangles()

        FloattingBody.__init__(self, combine.vertices, combine.faces)


class RectangularParallelepiped(FloattingBody):
    def __init__(self, height=10.0, length=10.0, thickness=2.0, nh=5, nl=5, nth=3):
        sides = OpenRectangularParallelepiped(height=height, length=length, thickness=thickness, nh=nh, nl=nl, nth=nth)
        top = OneSidedRectangle(height=thickness, length=length, nh=nth, nl=nl)
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


if __name__ == "__main__":
    from meshmagick.mmio import write_mesh

    rp = RectangularParallelepiped(height=10.0, length=10.0, thickness=10.0, nh=6, nl=6, nth=6)
    # rp.show_matplotlib()
    # write_mesh("test.dat", rp.vertices, rp.faces, "mar")
