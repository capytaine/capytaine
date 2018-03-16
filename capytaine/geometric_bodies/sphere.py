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


def generate_sphere(radius=1.0, ntheta=10, nphi=10,
                    z0=0.0, clip_free_surface=False, half=False,
                    name=None):
    """Generate the mesh of a sphere.

    Parameters
    ----------
    radius : float
        radius of the sphere
    ntheta : int
        number of panels along a meridian (or number of parallels-1)
    nphi : int
        number of panels along a parallel (or number of meridian-1)
    z0 : float
        depth of the center of mass of the sphere
    clip_free_surface : bool
        if True, only mesh the part of the sphere where z < 0,
        can be used with z0 to obtain any clipped sphere.
    half : bool
        if True, only mesh the part of the sphere where y > 0
    name : string
        a name identifying the sphere (default: "sphere_id" where id is an unique integer).

    Returns
    -------
    FloatingBody
        the generated body
    """

    if clip_free_surface:
        if z0 < -radius:  # fully immersed
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Sphere out of the water")
    else:
        theta_max = np.pi

    theta = np.linspace(0.0, theta_max, ntheta+1)
    if half:
        phi = np.linspace(-np.pi/2, np.pi/2, nphi+1)
    else:
        phi = np.linspace(-np.pi, np.pi, nphi+1)

    # Nodes
    nodes = np.zeros(((ntheta+1)*(nphi+1), 3), dtype=np.float32)

    for i, (t, p) in enumerate(product(theta, phi)):
        # The sign of theta below is a trick to get the correct orientation of the normal vectors...
        x = radius * np.sin(t) * np.sin(np.sign(t)*p)
        y = radius * np.sin(t) * np.cos(np.sign(t)*p)
        z = z0 - radius * np.cos(t)
        nodes[i, :] = (x, y, z)

    # Connectivity
    panels = np.zeros((ntheta*nphi, 4), dtype=np.int)

    for k, (i, j) in enumerate(product(range(0, ntheta), range(0, nphi))):
        panels[k, :] = (j+i*(nphi+1), j+(i+1)*(nphi+1), j+1+(i+1)*(nphi+1), j+1+i*(nphi+1))

    if name is None:
        name = f"sphere_{next(Mesh._ids)}"
    sphere = FloatingBody(nodes, panels, name=name)
    sphere.mesh.merge_duplicates()
    sphere.mesh.heal_triangles()

    return sphere


def generate_half_sphere(**kwargs):
    return generate_sphere(half=True, **kwargs)


def generate_clever_sphere(radius=1.0, ntheta=10, nphi=10,
                           z0=0.0, clip_free_surface=False,
                           name=None):
    """Generate the floating body of a sphere using its axial symmetry.

    Same arguments as `generate_sphere`."""
    if clip_free_surface:
        if z0 < -radius:  # fully immersed
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Sphere out of the water")
    else:
        theta_max = np.pi
    theta = np.linspace(0.0, theta_max, ntheta+1)
    circle_profile = np.empty((ntheta+1, 3), dtype=np.float32)
    for i, t in enumerate(theta):
        x = radius * np.sin(t)
        z = z0 - radius * np.cos(t)
        circle_profile[i, :] = (x, 0, z)

    if name is None:
        name = f"sphere_{next(Mesh._ids)}"
    return AxialSymmetry.from_profile(circle_profile, point_on_rotation_axis=np.zeros(3), nphi=nphi, name=name)

