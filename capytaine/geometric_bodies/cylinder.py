#!/usr/bin/env python
# coding: utf-8
"""
Generate meshes of cylinders and disks

This file is part of "capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
from itertools import product

import numpy as np

from meshmagick.mesh import Mesh
from meshmagick.geometry import xOz_Plane, yOz_Plane

from capytaine.bodies import FloatingBody
from capytaine.meshes_collection import CollectionOfMeshes
from capytaine.symmetries import TranslationalSymmetry, ReflectionSymmetry

LOG = logging.getLogger(__name__)


##########
#  Disk  #
##########

class Disk(FloatingBody):
    def __init__(self, radius=1.0, center=(0, 0, 0), nr=3, ntheta=5, name=None):
        self.radius = radius
        self.center = np.asarray(center)

        if name is None:
            name = f"disk_{next(Mesh._ids)}"

        mesh = self.generate_disk_mesh(nr, ntheta, name)

        FloatingBody.__init__(self, mesh=mesh, name=name)

    def generate_disk_mesh(self, nr: int, ntheta: int, name=None):
        theta_max = np.pi
        theta = np.linspace(-theta_max, theta_max, ntheta+1)
        R = np.linspace(0.0, self.radius, nr+1)

        nodes = np.zeros(((ntheta+1)*(nr+1), 3), dtype=np.float)

        for i, (r, t) in enumerate(product(R, theta)):
            y = +r * np.sin(t)
            z = -r * np.cos(t)
            nodes[i, :] = (0, y, z)
        nodes += self.center

        panels = np.zeros((ntheta*nr, 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, nr), range(0, ntheta))):
            panels[k, :] = (
                j+i*(ntheta+1),
                j+1+i*(ntheta+1),
                j+1+(i+1)*(ntheta+1),
                j+(i+1)*(ntheta+1)
            )

        mesh = Mesh(nodes, panels, name=f"{name}_mesh")
        mesh.merge_duplicates()
        mesh.heal_triangles()

        return mesh


##############
#  Cylinder  #
##############

class HorizontalCylinder(FloatingBody):
    def __init__(self, length=10.0, radius=1.0, center=(0, 0, 0),
                 nx=10, ntheta=10, nr=2,
                 clever=True, name=None):
        """Generate the mesh of an horizontal cylinder.

        Parameters
        ----------
        length : float
            length of the cylinder
        radius : float
            radius of the cylinder
        center : 3-ple or array of shape (3,)
            position of the center of the cylinder
        nx : int
            number of circular slices
        ntheta : int
            number of panels along a circular slice of the cylinder
        clip_free_surface : bool
            if True, only mesh the part of the cylinder where z < 0,
            can be used with z0 to obtain any clipped cylinder
        """
        self.length = length
        self.radius = radius
        self.center = np.asarray(center, dtype=np.float)

        if name is None:
            name = f"cylinder_{next(Mesh._ids)}"

        open_cylinder = self._generate_open_cylinder_mesh(nx, ntheta, name)

        if nr > 0:
            side = Disk(radius=radius, center=(-np.array([length/2, 0, 0])),
                        nr=nr, ntheta=ntheta, name=f"side_of_{name}").mesh

            other_side = side.copy()
            other_side.name = f"other_side_of_{name}"
            other_side.mirror(yOz_Plane)

            mesh = CollectionOfMeshes((open_cylinder, side, other_side))

        else:
            mesh = open_cylinder

        if not clever:
            mesh = mesh.merge()
            mesh.merge_duplicates()
            mesh.heal_triangles()

        mesh.translate(self.center)

        FloatingBody.__init__(self, mesh=mesh, name=name)

    def _generate_open_cylinder_mesh(self, nx, ntheta, name=None):
        """Open horizontal cylinder using the symmetry to speed up the computations"""
        theta_max = np.pi
        theta = np.linspace(0, theta_max, ntheta//2+1)
        X = np.array([0, self.length/nx])

        # Nodes
        nodes = np.zeros(((ntheta//2+1)*2, 3), dtype=np.float)

        for i, (t, x) in enumerate(product(theta, X)):
            y = + self.radius * np.sin(t)
            z = - self.radius * np.cos(t)
            nodes[i, :] = (x, y, z)
        nodes += -np.array([self.length/2, 0, 0])

        # Connectivities
        panels = np.zeros((ntheta//2, 4), dtype=np.int)

        for k, i in enumerate(range(0, ntheta//2)):
            panels[k, :] = (2*i, 2*i+2, 2*i+3, 2*i+1)
        half_ring = Mesh(nodes, panels, name=f"half_ring_of_{name}_mesh")

        ring = ReflectionSymmetry(half_ring, plane=xOz_Plane, name=f"ring_of_{name}_mesh")

        return TranslationalSymmetry(ring, translation=np.asarray([self.length/nx, 0.0, 0.0]),
                                     nb_repetitions=nx-1, name=f"{name}_mesh")

    @property
    def volume(self):
        return self.length*np.pi*self.radius**2

