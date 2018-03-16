#!/usr/bin/env python
# coding: utf-8
"""
Generate mesh for rectangles and parallelepipeds.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

from itertools import product

import numpy as np

from meshmagick.mesh import Mesh
from capytaine.bodies import FloatingBody
from capytaine.symmetries import TranslationalSymmetry


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


