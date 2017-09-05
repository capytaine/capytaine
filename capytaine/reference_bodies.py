#!/usr/bin/env python
# coding: utf-8
"""
class Sphere
class HorizontalCylinder
class OneSidedRectangle
class TwoSidedRectangle
class OpenRectangularParallelepiped
class RectangularParallelepiped
"""

from itertools import product
import numpy as np

from capytaine.bodies import *

class Sphere(FloattingBody):
    """Floatting body of the shape of a sphere."""

    def __init__(self,
            radius=1.0, ntheta=11, nphi=11,
            z0=0.0, clip_free_surface=False,
            half=False):

        if clip_free_surface:
            if z0 < -radius: # fully immerged
                theta_max = np.pi
            elif z0 < radius:
                theta_max = np.arccos(z0/radius)
            else:
                raise Exception("Sphere out of the water")
        else:
            theta_max = np.pi

        if half:
            theta = np.linspace(-theta_max, 0.0, ntheta)
        else:
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
        self.merge_duplicates()
        self.heal_triangles()


class HalfSphere(FloattingBody):
    """Floatting body of the shape of half a sphere."""

    def __init__(self, **kwargs):
        Sphere.__init__(self, half=True, **kwargs)


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

    def __init__(self, height=2.0, width=10.0, nh=5, nw=5, z0=0.0):

        X = np.linspace(-width/2, width/2, nw)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nw*nh, 3), dtype=np.float32)
        panels = np.zeros(((nw-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)

        FloattingBody.__init__(self, nodes, panels)


class TwoSidedRectangle(FloattingBody):
    """Rectangular panel with cartesian mesh.
    Each face is defined twice with two opposite normal vectors."""

    def __init__(self, height=2.0, width=10.0, nh=5, nw=5, z0=0.0):

        X = np.linspace(-width/2, width/2, nw)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nw*nh, 3), dtype=np.float32)
        panels = np.zeros((2*(nw-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)
            panels[(nw-1)*(nh-1)+k, :] = (j+i*nh, j+(i+1)*nh, j+1+(i+1)*nh, j+1+i*nh)

        FloattingBody.__init__(self, nodes, panels)


class OpenRectangularParallelepiped(FloattingBody):
    """Four panels forming a parallelepiped without top nor bottom."""

    def __init__(self, height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=3, z0=0.0):
        front = OneSidedRectangle(height=height, width=width, nh=nh, nw=nw, z0=z0)
        back = front.copy()

        front.translate_y(thickness/2)
        back.rotate_z(np.pi)
        back.translate_y(-thickness/2)

        side = OneSidedRectangle(height=height, width=thickness, nh=nh, nw=nth, z0=z0)
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

    def __init__(self, height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=3):
        sides = OpenRectangularParallelepiped(height=height, width=width, thickness=thickness, nh=nh, nw=nw, nth=nth)
        top = OneSidedRectangle(height=thickness, width=width, nh=nth, nw=nw)
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
