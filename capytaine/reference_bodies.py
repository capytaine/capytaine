#!/usr/bin/env python
# coding: utf-8
"""
Generate mesh for some simple geometric shapes.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

from itertools import product

import numpy as np

from meshmagick.mesh import Mesh
from capytaine.bodies import FloatingBody
from capytaine.symmetries import xOz_Plane, yOz_Plane, ReflectionSymmetry, TranslationalSymmetry, AxialSymmetry


##########
#  Disk  #
##########

def generate_disk(radius=1.0, nr=3, ntheta=5,
                  z0=0.0, clip_free_surface=False,
                  name=None):

    if clip_free_surface:
        if z0 < -radius: # fully immersed
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
    disk = FloatingBody(nodes, panels, name=name)
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
        if z0 < -radius: # fully immersed
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
    cylinder = FloatingBody(nodes, panels, name=name)
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


################
#  Rectangles  #
################

def generate_one_sided_rectangle(height=5.0, width=5.0, nh=5, nw=5, name=None):
    """Generate the mesh of a rectangle.

    Normals are oriented in the positive y direction.

    Parameters
    ----------
    height : float
        height of the panel (size along z)
    width : float
        width of the panel (size along x)
    nh : int
        number of panels in the z direction
    nw : int
        number of panels in the x direction

    Returns
    -------
    FloatingBody
        the generated body
    """

    X = np.linspace(-width/2, width/2, nw+1)
    Z = np.linspace(0, height, nh+1)

    nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=np.float32)
    panels = np.zeros((nw*nh, 4), dtype=np.int)

    for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
        nodes[i, :] = x, y, z

    for k, (i, j) in enumerate(product(range(0, nw), range(0, nh))):
        panels[k, :] = (j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))

    if name is None:
        name = f"rectangle_{next(Mesh._ids)}"
    return FloatingBody(nodes, panels, name=name)


def generate_clever_one_sided_rectangle(width=5.0, nw=5, name=None, **kwargs):
    if name is None:
        name = f"rectangle_{next(Mesh._ids)}"
    strip = generate_one_sided_rectangle(width=width/nw, nw=1, name=f"strip_of_{name}", **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1, name=name)


def generate_free_surface(width=100, length=100, nw=10, nl=10, name=None):
    """Special body for the meshing of the free surface."""
    X = np.linspace(-width/2, width/2, nw+1)
    Y = np.linspace(-length/2, length/2, nl+1)

    nodes = np.zeros(((nw+1)*(nl+1), 3), dtype=np.float32)
    panels = np.zeros((nw*nl, 4), dtype=np.int)

    for i, (x, y, z) in enumerate(product(X, Y, [0.0])):
        nodes[i, :] = x, y, z

    for k, (i, j) in enumerate(product(range(0, nw), range(0, nl))):
        panels[k, :] = (j+i*(nl+1), j+1+i*(nl+1), j+1+(i+1)*(nl+1), j+(i+1)*(nl+1))

    if name is None:
        name = f"free_surface_{next(Mesh._ids)}"
    return FloatingBody(nodes, panels, name=name)


def generate_open_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0,
                                             nh=5, nw=5, nth=1,
                                             name=None):
    """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

    Parameters
    ----------
    height : float
        height of the object (size along z)
    width : float
        width of the object (size along x)
    thickness : float
        thickness of the object (size along y)
    nh : int
        number of panels in the z direction
    nw : int
        number of panels in the x direction
    nth : int
        number of panels in the y direction

    Returns
    -------
    FloatingBody
        the generated body
    """

    if name is None:
        name = f"open_parallelepiped_{next(Mesh._ids)}"

    front = generate_one_sided_rectangle(height=height, width=width, nh=nh, nw=nw, name=f"front_of_{name}")
    back = front.copy(name=f"back_of_{name}")

    front.translate_y(thickness/2)
    back.rotate_z(np.pi)
    back.translate_y(-thickness/2)

    parallelepiped = front + back

    if nth > 0:
        side = generate_one_sided_rectangle(height=height, width=thickness, nh=nh, nw=nth, name=f"side_of_{name}")
        other_side = side.copy(name=f"other_side_of_{name}")

        side.rotate_z(np.pi/2)
        side.translate_x(-width/2)
        other_side.rotate_z(-np.pi/2)
        other_side.translate_x(width/2)

        parallelepiped = parallelepiped + side + other_side

    parallelepiped = parallelepiped.as_FloatingBody(name=name)
    parallelepiped.mesh.merge_duplicates()
    parallelepiped.mesh.heal_triangles()

    return parallelepiped


def generate_clever_open_rectangular_parallelepiped(width=5.0, nw=5, name=None, **kwargs):
    if name is None:
        name = f"open_parallelepiped_{next(Mesh._ids)}"
    strip = generate_open_rectangular_parallelepiped(width=width/nw, nw=1, nth=0, name=f"strip_of_{name}", **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1, name=name)


def generate_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=1, name=None):
    """Generate the mesh of six rectangles forming a complete rectangular parallelepiped.

    Parameters
    ----------
    height : float
        height of the object (size along z)
    width : float
        width of the object (size along x)
    thickness : float
        thickness of the object (size along y)
    nh : int
        number of panels in the z direction
    nw : int
        number of panels in the x direction
    nth : int
        number of panels in the y direction
    z0 : float
        depth of the bottom of the object

    Returns
    -------
    FloatingBody
        the generated body
    """

    if name is None:
        name = f"parallelepiped_{next(Mesh._ids)}"

    sides = generate_open_rectangular_parallelepiped(
        height=height, width=width, thickness=thickness,
        nh=nh, nw=nw, nth=nth,
        name=f"sides_of_{name}")

    top = generate_one_sided_rectangle(
        height=thickness, width=width,
        nh=nth, nw=nw,
        name=f"top_of_{name}")
    bottom = top.copy(name=f"bottom_of_{name}")

    top.rotate_x(np.pi/2)
    top.translate_y(thickness/2)
    top.translate_z(height)
    bottom.rotate_x(-np.pi/2)
    bottom.translate_y(-thickness/2)

    parallelepiped = sides + top + bottom
    parallelepiped = parallelepiped.as_FloatingBody(name=name)
    parallelepiped.mesh.merge_duplicates()
    parallelepiped.mesh.heal_triangles()

    return parallelepiped


def generate_horizontal_open_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0,
                                                        nh=5, nw=5, nth=1,
                                                        **kwargs):
    orp = generate_open_rectangular_parallelepiped(
        height=width, width=height, thickness=thickness,
        nh=nw, nw=nh, nth=nth,
        **kwargs)
    orp.rotate_y(-np.pi/2)
    return orp


def generate_clever_horizontal_open_rectangular_parallelepiped(width=10.0, nw=5, name=None, **kwargs):
    if name is None:
        name = f"open_parallelepiped_{next(Mesh._ids)}"
    strip = generate_horizontal_open_rectangular_parallelepiped(width=width/nw, nw=1, name=f"strip_of_{name}", **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1, name=name)


###########
#  Other  #
###########

def generate_dummy_floating_body():
    return FloatingBody(np.zeros((0, 3)), np.zeros((0, 4)), name="dummy_body")
