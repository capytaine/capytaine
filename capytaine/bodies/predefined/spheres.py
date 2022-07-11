#!/usr/bin/env python
# coding: utf-8
"""Generate spherical bodies."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np

from capytaine.meshes.geometry import Axis
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric import AxialSymmetricMesh
from capytaine.bodies.bodies import FloatingBody

LOG = logging.getLogger(__name__)


class Sphere(FloatingBody):
    """Sphere

    Parameters
    ----------
    radius : float
        radius of the sphere
    center : 3-ple or array of shape (3,)
        position of the geometric center of the sphere
    ntheta : int
        number of panels along a meridian (or number of parallels-1)
    nphi : int
        number of panels along a parallel (or number of meridians-1)
    axial_symmetry : bool
        if True, use the axial symmetry to build the mesh (default: True)
    clip_free_surface : bool
        if True, only mesh the part of the sphere where z < 0 (default: False),
        can be used with center to obtain any clipped sphere,
        if True, then ntheta is the number of parallel below the free surface.
    name : string
        a name identifying the sphere (default: "sphere_id" where id is an unique integer).
    """

    def __init__(self, *, radius=1.0, center=(0, 0, 0),
                 ntheta=10, nphi=10, clip_free_surface=False,
                 axial_symmetry=True, clever=None,
                 name=None):
        self.radius = radius
        self.geometric_center = np.array(center, dtype=float)

        if clever is not None:
            LOG.warning("Deprecation warning: `clever` argument for Sphere is deprecated. "
                        "Use `axial_symmetry` instead.")

        if name is None:
            name = f"sphere_{next(Mesh._ids)}"

        mesh = self._generate_mesh_using_symmetry(ntheta, nphi, clip_free_surface, f"{name}_mesh")

        if not axial_symmetry:
            mesh = mesh.merged()
            mesh.merge_duplicates()
            mesh.heal_triangles()

        FloatingBody.__init__(self, mesh=mesh, name=name)


    def _generate_mesh_using_symmetry(self, ntheta, nphi, clip_free_surface, mesh_name):
        if clip_free_surface:
            if self.geometric_center[2] < -self.radius:  # fully immersed
                theta_max = np.pi
            elif self.geometric_center[2] < self.radius:
                theta_max = np.arccos(self.geometric_center[2]/self.radius)
            else:
                raise ValueError("Impossible to mesh the immersed hull of a sphere completely out of the water")
        else:
            theta_max = np.pi

        theta = np.linspace(0.0, theta_max, ntheta+1)
        points_on_a_meridian = (
                self.radius * np.stack([np.sin(theta), np.zeros_like(theta), -np.cos(theta)], axis=1)
                + self.geometric_center
                )

        symmetry_axis = Axis(vector=[0, 0, 1], point=self.geometric_center)
        return AxialSymmetricMesh.from_profile(points_on_a_meridian, axis=symmetry_axis, nphi=nphi, name=mesh_name)

