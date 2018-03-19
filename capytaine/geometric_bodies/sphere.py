#!/usr/bin/env python
# coding: utf-8
"""
Generate meshes of spheres

This file is part of "capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
from itertools import product

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.bodies import FloatingBody
from capytaine.symmetries import AxialSymmetry

LOG = logging.getLogger(__name__)


class Sphere(FloatingBody):
    def __init__(self, radius=1.0, center=(0, 0, 0),
                 ntheta=10, nphi=10, clever=True, clip_free_surface=False,
                 name=None):
        """Generate a sphere.

        Parameters
        ----------
        radius : float
            radius of the sphere
        center : 3-ple or array of shape (3,)
            position of the center of sphere
        ntheta : int
            number of panels along a meridian (or number of parallels-1)
        nphi : int
            number of panels along a parallel (or number of meridian-1)
        clever : bool
            if True, use the symmetries to build the mesh (default: True)
        clip_free_surface : bool
            if True, only mesh the part of the sphere where z < 0 (default: False),
            can be used with center to obtain any clipped sphere.
        name : string
            a name identifying the sphere (default: "sphere_id" where id is an unique integer).
        """
        self.radius = radius
        self.center = np.asarray(center)

        if not clever:
            mesh = self.generate_sphere_mesh(ntheta, nphi, clip_free_surface, name)
        else:
            mesh = self.generate_clever_sphere_mesh(ntheta, nphi, clip_free_surface, name)

        if name is None:
            name = f"sphere_{next(Mesh._ids)}"

        FloatingBody.__init__(self, mesh, name)

    def generate_sphere_mesh(self, ntheta=10, nphi=10, clip_free_surface=False, name=None):
        if clip_free_surface:
            if self.center[2] < -self.radius:  # fully immersed
                theta_max = np.pi
            elif self.center[2] < self.radius:
                theta_max = np.arccos(self.center[2]/self.radius)
            else:
                raise Exception("Sphere out of the water")
        else:
            theta_max = np.pi

        theta = np.linspace(0.0, theta_max, ntheta+1)
        phi = np.linspace(-np.pi, np.pi, nphi+1)

        # Nodes
        nodes = np.zeros(((ntheta+1)*(nphi+1), 3), dtype=np.float32)

        for i, (t, p) in enumerate(product(theta, phi)):
            # The sign of theta below is a trick to get the correct orientation of the normal vectors...
            x = + np.sin(t) * np.sin(np.sign(t)*p)
            y = + np.sin(t) * np.cos(np.sign(t)*p)
            z = - np.cos(t)
            nodes[i, :] = (x, y, z)
        nodes *= self.radius
        nodes += self.center

        # Connectivity
        panels = np.zeros((ntheta*nphi, 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta), range(0, nphi))):
            panels[k, :] = (j+i*(nphi+1), j+(i+1)*(nphi+1), j+1+(i+1)*(nphi+1), j+1+i*(nphi+1))

        mesh = Mesh(nodes, panels, name=f"{name}_mesh")
        mesh.merge_duplicates()
        mesh.heal_triangles()

        return mesh

    def generate_clever_sphere_mesh(self, ntheta=10, nphi=10, clip_free_surface=False, name=None):
        if clip_free_surface:
            if self.center[2] < -self.radius:  # fully immersed
                theta_max = np.pi
            elif self.center[2] < self.radius:
                theta_max = np.arccos(self.center[2]/self.radius)
            else:
                raise Exception("Sphere out of the water")
        else:
            theta_max = np.pi

        theta = np.linspace(0.0, theta_max, ntheta+1)

        circle_profile = np.zeros((ntheta+1, 3), dtype=np.float32)
        circle_profile[:, 0] = np.sin(theta)
        circle_profile[:, 2] = -np.cos(theta)
        circle_profile *= self.radius
        circle_profile += self.center

        return AxialSymmetry.from_profile(circle_profile, point_on_rotation_axis=self.center,
                                          nphi=nphi, name=f"{name}_mesh")

    @property
    def volume(self):
        return 4/3*np.pi*self.radius**3
