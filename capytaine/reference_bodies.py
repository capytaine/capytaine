#!/usr/bin/env python
# coding: utf-8
"""
Generate mesh for some simple geometric shapes.
"""

from itertools import product, count

import numpy as np

from capytaine.bodies import FloatingBody
from capytaine.symmetries import TranslationalSymmetry, AxialSymmetry

#############
#  Spheres  #
#############

def generate_sphere(radius=1.0, ntheta=10, nphi=10,
                    z0=0.0, clip_free_surface=False, half=False):
    """Generate the mesh of a sphere.

    Parameters
    ----------
    radius: float
        radius of the sphere
    ntheta: int
        number of panels along a meridian (or number of parallels-1)
    nphi: int
        number of panels along a parallel (or number of meridian-1)
    z0: float
        depth of the center of mass of the sphere
    clip_free_surface: bool
        if True, only mesh the part of the sphere where z < 0,
        can be used with z0 to obtain any clipped sphere.
    half: bool
        if True, only mesh the part of the sphere where y > 0
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

    if half:
        theta = np.linspace(0.0, theta_max, ntheta+1)
    else:
        theta = np.linspace(-theta_max, theta_max, ntheta+1)
    phi = np.linspace(-np.pi/2, np.pi/2, nphi+1)

    # Nodes
    nodes = np.zeros(((ntheta+1)*(nphi+1), 3), dtype=np.float32)

    for i, (t, p) in enumerate(product(theta, phi)):
        # The sign of theta below is a trick to get the correct orientation of the normal vectors...
        x =      radius * np.sin(t) * np.sin(np.sign(t)*p)
        y =      radius * np.sin(t) * np.cos(np.sign(t)*p)
        z = z0 - radius * np.cos(t)
        nodes[i, :] = (x, y, z)

    # Connectivities
    panels = np.zeros((ntheta*nphi, 4), dtype=np.int)

    for k, (i, j) in enumerate(product(range(0, ntheta), range(0, nphi))):
        panels[k, :] = (j+i*(nphi+1), j+(i+1)*(nphi+1), j+1+(i+1)*(nphi+1), j+1+i*(nphi+1))

    sphere = FloatingBody(nodes, panels, name=f"sphere_{next(FloatingBody._ids)}")
    sphere.merge_duplicates()
    sphere.heal_triangles()

    return sphere


def generate_half_sphere(**kwargs):
    return generate_sphere(half=True, **kwargs)


def generate_axi_symmetric_body(profile, point_on_rotation_axis=np.zeros(3), nth=20):
    profile = np.asarray(profile)
    assert len(profile.shape) == 2
    assert profile.shape[1] == 3

    n = profile.shape[0]
    angle = 2*np.pi/nth

    rotated_profile = FloatingBody(profile, np.zeros((0, 4)))
    rotated_profile.rotate_z(angle)

    nodes_slice = np.concatenate([profile, rotated_profile.vertices])
    faces_slice = np.array([[i, i+n, i+n+1, i+1] for i in range(n-1)])
    body_slice = FloatingBody(nodes_slice, faces_slice)
    body_slice.merge_duplicates()
    body_slice.heal_triangles()

    return AxialSymmetry(body_slice, point_on_rotation_axis=point_on_rotation_axis, nb_repetitions=nth-1)


###############
#  Cylinders  #
###############


def generate_horizontal_cylinder(length=10.0, radius=1.0,
                                 nx=10, nr=2, ntheta=10,
                                 z0=0.0, clip_free_surface=False):
    """Generate the mesh of an horizontal cylinder.

    Parameters
    ----------
    length: float
        length of the cylinder
    radius: float
        radius of the cylinder
    nx: int
        number of circular slices
    nr: int
        at the ends of the cylinder, number of panels along a radius
    ntheta: int
        number of panels along a circular slice of the cylinder
    z0: float
        depth of the bottom of the cylinder
    clip_free_surface: bool
        if True, only mesh the part of the cylinder where z < 0,
        can be used with z0 to obtain any clipped cylinder
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

    theta = np.linspace(-theta_max, theta_max, ntheta+1)
    X = np.linspace(0.0, length, nx+1)
    R = np.linspace(0.0, radius, nr+1)

    # Nodes
    nodes = np.zeros(((ntheta+1)*(nx+2*nr+3), 3), dtype=np.float32)

    for i, (t, x) in enumerate(product(theta, X)):
        y = radius * np.sin(t)
        z = z0 - radius * np.cos(t)
        nodes[i, :] = (x, y, z)

    for i, (x, r, t) in enumerate(product([0, length], R, theta)):
        y = r * np.sin(t)
        z = z0 - r * np.cos(t)
        nodes[(ntheta+1)*(nx+1)+i, :] = (x, y, z)

    # Connectivities
    panels = np.zeros((ntheta*(nx+2*nr), 4), dtype=np.int)

    for k, (i, j) in enumerate(product(range(0, ntheta),
                                       range(0, nx))):
        panels[k, :] = (
            j+i*(nx+1),
            j+(i+1)*(nx+1),
            j+1+(i+1)*(nx+1),
            j+1+i*(nx+1)
        )

    for k, (i, j) in enumerate(product(range(0, nr),
                                       range((ntheta+1)*(nx+1), (ntheta+1)*(nx+1)+ntheta))):
        panels[ntheta*nx+k, :] = (
            j+i*(ntheta+1),
            j+1+i*(ntheta+1),
            j+1+(i+1)*(ntheta+1),
            j+(i+1)*(ntheta+1)
        )

    for k, (i, j) in enumerate(product(range(0, nr),
                                       range((ntheta+1)*((nx+1)+(nr+1)), (ntheta+1)*((nx+1)+(nr+1))+ntheta))):
        panels[ntheta*(nx+nr)+k, :] = (
            j+i*(ntheta+1),
            j+(i+1)*(ntheta+1),
            j+1+(i+1)*(ntheta+1),
            j+1+i*(ntheta+1)
        )

    cylinder = FloatingBody(nodes, panels, name=f"cylinder_{next(FloatingBody._ids)}")
    cylinder.merge_duplicates()
    cylinder.heal_triangles()

    return cylinder

def generate_ring(**kwargs):
    return generate_horizontal_cylinder(nx=1, **kwargs)

def generate_clever_horizontal_cylinder(length=10, nx=10, **kwargs):
    """Open horizontal cylinder using the symmetry to speed up the computations"""
    ring = generate_ring(length=length/nx, nr=0, **kwargs)
    return TranslationalSymmetry(ring, translation=np.asarray([length/nx, 0.0, 0.0]), nb_repetitions=nx-1)


################
#  Rectangles  #
################

def generate_one_sided_rectangle(height=5.0, width=5.0, nh=5, nw=5):
    """Generate the mesh of a rectangle.

    Normals are oriented in the positive y direction.

    Parameters
    ----------
    height: float
        height of the panel (size along z)
    width: float
        width of the panel (size along x)
    nh: int
        number of panels in the z direction
    nw: int
        number of panels in the x direction
    """

    X = np.linspace(-width/2, width/2, nw+1)
    Z = np.linspace(0, height, nh+1)

    nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=np.float32)
    panels = np.zeros((nw*nh, 4), dtype=np.int)

    for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
        nodes[i, :] = x, y, z

    for k, (i, j) in enumerate(product(range(0, nw), range(0, nh))):
        panels[k, :] = (j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))

    return FloatingBody(nodes, panels, name=f"rectangle_{next(FloatingBody._ids)}")

def generate_clever_one_sided_rectangle(width=5.0, nw=5, **kwargs):
    strip = generate_one_sided_rectangle(width=width/nw, nw=1, **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1)

def generate_free_surface(width=100, length=100, nw=10, nl=10):
    """ """
    X = np.linspace(-width/2, width/2, nw+1)
    Y = np.linspace(-length/2, length/2, nl+1)

    nodes = np.zeros(((nw+1)*(nl+1), 3), dtype=np.float32)
    panels = np.zeros((nw*nl, 4), dtype=np.int)

    for i, (x, y, z) in enumerate(product(X, Y, [0.0])):
        nodes[i, :] = x, y, z

    for k, (i, j) in enumerate(product(range(0, nw), range(0, nl))):
        panels[k, :] = (j+i*(nl+1), j+1+i*(nl+1), j+1+(i+1)*(nl+1), j+(i+1)*(nl+1))

    return FloatingBody(nodes, panels, name=f"free_surface_{next(FloatingBody._ids)}")

def generate_open_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0,
                                             nh=5, nw=5, nth=1):
    """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

    Parameters
    ----------
    height: float
        height of the object (size along z)
    width: float
        width of the object (size along x)
    thickness: float
        thickness of the object (size along y)
    nh: int
        number of panels in the z direction
    nw: int
        number of panels in the x direction
    nth: int
        number of panels in the y direction
    """
    front = generate_one_sided_rectangle(height=height, width=width, nh=nh, nw=nw)
    back = front.copy()

    front.translate_y(thickness/2)
    back.rotate_z(np.pi)
    back.translate_y(-thickness/2)

    parallelepiped = front + back

    if nth > 0:
        side = generate_one_sided_rectangle(height=height, width=thickness, nh=nh, nw=nth)
        other_side = side.copy()

        side.rotate_z(np.pi/2)
        side.translate_x(-width/2)
        other_side.rotate_z(-np.pi/2)
        other_side.translate_x(width/2)

        parallelepiped = parallelepiped + side + other_side

    parallelepiped = parallelepiped.as_FloatingBody()
    parallelepiped.name = f"parallelepiped_{next(FloatingBody._ids)}"
    parallelepiped.merge_duplicates()
    parallelepiped.heal_triangles()

    return parallelepiped

def generate_clever_open_rectangular_parallelepiped(width=5.0, nw=5, **kwargs):
    strip = generate_open_rectangular_parallelepiped(width=width/nw, nw=1, nth=0, **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1)

def generate_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=1):
    """Generate the mesh of six rectangles forming a complete rectangular parallelepiped.
    Parameters
    ----------
    height: float
        height of the object (size along z)
    width: float
        width of the object (size along x)
    thickness: float
        thickness of the object (size along y)
    nh: int
        number of panels in the z direction
    nw: int
        number of panels in the x direction
    nth: int
        number of panels in the y direction
    z0: float
        depth of the bottom of the object
    """
    sides = generate_open_rectangular_parallelepiped(
        height=height, width=width, thickness=thickness,
        nh=nh, nw=nw, nth=nth)

    top = generate_one_sided_rectangle(
        height=thickness, width=width,
        nh=nth, nw=nw)
    bottom = top.copy()

    top.rotate_x(np.pi/2)
    top.translate_y(thickness/2)
    top.translate_z(height)
    bottom.rotate_x(-np.pi/2)
    bottom.translate_y(-thickness/2)

    parallelepiped = sides + top + bottom
    parallelepiped = parallelepiped.as_FloatingBody()
    parallelepiped.name = f"parallelepiped_{next(FloatingBody._ids)}"
    parallelepiped.merge_duplicates()
    parallelepiped.heal_triangles()

    return parallelepiped

def generate_horizontal_open_rectangular_parallelepiped(height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=1):
    orp = generate_open_rectangular_parallelepiped(
        height=width, width=height, thickness=thickness,
        nh=nw, nw=nh, nth=nth)
    orp.rotate_y(-np.pi/2)
    return orp

def generate_clever_horizontal_open_rectangular_parallelepiped(width=10.0, nw=5, **kwargs):
    strip = generate_horizontal_open_rectangular_parallelepiped(width=width/nw, nw=1, **kwargs)
    return TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1)

###########
#  Other  #
###########

def generate_dummy_floating_body():
    return FloatingBody(np.zeros((0, 3)), np.zeros((0, 4)), name="dummy_body")
