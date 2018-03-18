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
from capytaine.bodies import FloatingBody
from capytaine.symmetries import xOz_Plane, yOz_Plane, TranslationalSymmetry, ReflectionSymmetry

LOG = logging.getLogger(__name__)


##########
#  Disk  #
##########

def generate_disk(radius=1.0, nr=3, ntheta=5,
                  z0=0.0, clip_free_surface=False,
                  name=None):

    if clip_free_surface:
        if z0 < -radius:  # fully immersed
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Disk out of the water")
    else:
        theta_max = np.pi

    theta = np.linspace(-theta_max, theta_max, ntheta+1)
    R = np.linspace(0.0, radius, nr+1)

    nodes = np.zeros(((ntheta+1)*(nr+1), 3), dtype=np.float32)

    for i, (r, t) in enumerate(product(R, theta)):
        y = r * np.sin(t)
        z = z0 - r * np.cos(t)
        nodes[i, :] = (0, y, z)

    panels = np.zeros((ntheta*nr, 4), dtype=np.int)

    for k, (i, j) in enumerate(product(range(0, nr), range(0, ntheta))):
        panels[k, :] = (
            j+i*(ntheta+1),
            j+1+i*(ntheta+1),
            j+1+(i+1)*(ntheta+1),
            j+(i+1)*(ntheta+1)
        )

    if name is None:
        name = f"disk_{next(Mesh._ids)}"
    disk = FloatingBody(Mesh(nodes, panels, name=f"{name}_mesh"), name=name)
    disk.mesh.merge_duplicates()
    disk.mesh.heal_triangles()

    return disk


###############
#  Cylinders  #
###############

def generate_open_horizontal_cylinder(length=10.0, radius=1.0,
                                      nx=10, ntheta=10,
                                      z0=0.0, clip_free_surface=False,
                                      half=False,
                                      name=None):
    """Generate the mesh of an horizontal cylinder.

    Parameters
    ----------
    length : float
        length of the cylinder
    radius : float
        radius of the cylinder
    nx : int
        number of circular slices
    ntheta : int
        number of panels along a circular slice of the cylinder
    z0 : float
        depth of the bottom of the cylinder
    clip_free_surface : bool
        if True, only mesh the part of the cylinder where z < 0,
        can be used with z0 to obtain any clipped cylinder
    half : bool
        if True, only mesh the part of the cylinder where y > 0

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
            raise Exception("Cylinder out of the water")
    else:
        theta_max = np.pi

    if half:
        theta = np.linspace(0.0, theta_max, ntheta+1)
    else:
        theta = np.linspace(-theta_max, theta_max, ntheta+1)
    X = np.linspace(0.0, length, nx+1)

    # Nodes
    nodes = np.zeros(((ntheta+1)*(nx+1), 3), dtype=np.float32)

    for i, (t, x) in enumerate(product(theta, X)):
        y = radius * np.sin(t)
        z = z0 - radius * np.cos(t)
        nodes[i, :] = (x, y, z)

    # Connectivities
    panels = np.zeros((ntheta*nx, 4), dtype=np.int)

    for k, (i, j) in enumerate(product(range(0, ntheta),
                                       range(0, nx))):
        panels[k, :] = (
            j+i*(nx+1),
            j+(i+1)*(nx+1),
            j+1+(i+1)*(nx+1),
            j+1+i*(nx+1)
        )

    if name is None:
        name = f"cylinder_{next(Mesh._ids)}"
    cylinder = FloatingBody(Mesh(nodes, panels, name=f"{name}_mesh"), name=name)
    cylinder.mesh.merge_duplicates()
    cylinder.mesh.heal_triangles()

    return cylinder


def generate_ring(**kwargs):
    if 'name' not in kwargs:
        kwargs['name'] = f"ring_{next(Mesh._ids)}"
    return generate_open_horizontal_cylinder(nx=1, **kwargs)


def generate_clever_horizontal_cylinder(length=10, nx=10, name=None, ntheta=10, **kwargs):
    """Open horizontal cylinder using the symmetry to speed up the computations"""
    if name is None:
        name = f"horizontal_cylinder_{next(Mesh._ids)}"
    half_ring = generate_ring(length=length/nx, name=f"half_slice_of_{name}", half=True, ntheta=ntheta//2, **kwargs)
    ring = ReflectionSymmetry(half_ring, plane=xOz_Plane)
    return TranslationalSymmetry(ring, translation=np.asarray([length/nx, 0.0, 0.0]), nb_repetitions=nx-1, name=name)


def generate_horizontal_cylinder(length=10.0, radius=1.0,
                                 nx=10, nr=2, ntheta=10,
                                 z0=0.0, clip_free_surface=False,
                                 name=None):
    """Generate the mesh of a closed horizontal cylinder.

    Parameters
    ----------
    length : float
        length of the cylinder
    radius : float
        radius of the cylinder
    nx : int
        number of circular slices
    nr : int
        at the ends of the cylinder, number of panels along a radius
    ntheta : int
        number of panels along a circular slice of the cylinder
    z0 : float
        depth of the bottom of the cylinder
    clip_free_surface : bool
        if True, only mesh the part of the cylinder where z < 0,
        can be used with z0 to obtain any clipped cylinder

    Returns
    -------
    FloatingBody
        the generated body
    """

    if name is None:
        name = f"cylinder_{next(Mesh._ids)}"

    open_cylinder = generate_open_horizontal_cylinder(
        length=length, radius=radius,
        nx=nx, ntheta=ntheta,
        z0=z0, clip_free_surface=clip_free_surface,
        name=f"body_of_{name}"
    )

    side = generate_disk(
        radius=radius,
        nr=nr, ntheta=ntheta,
        z0=z0, clip_free_surface=clip_free_surface,
        name=f"side_of_{name}"
    )

    other_side = side.copy(name=f"other_side_of_{name}")
    other_side.mirror(yOz_Plane)
    other_side.translate_x(length)

    cylinder = open_cylinder + side + other_side

    cylinder = cylinder.as_FloatingBody(name=name)
    cylinder.mesh.merge_duplicates()
    cylinder.mesh.heal_triangles()

    return cylinder
