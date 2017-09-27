#!/usr/bin/env python
# coding: utf-8
"""
class Sphere
class HorizontalCylinder
class OneSidedRectangle
class TwoSidedRectangle
class OpenRectangularParallelepiped
class RectangularParallelepiped
class DummyBody

TODO: Do we really need class and not just generating functions?
"""

from itertools import product, count

import numpy as np

from capytaine.bodies import FloatingBody


class Sphere(FloatingBody):
    """Floatting body of the shape of a sphere."""

    _ids = count(0)

    def __init__(self,
                 radius=1.0, ntheta=11, nphi=11,
                 z0=0.0, clip_free_surface=False,
                 half=False
                ):
        """Generate the mesh.

        Parameters
        ----------
        radius: float
            radius of the sphere
        ntheta: int
            number of points along a meridian (or number of parallels)
        nphi: int
            number of points along a parallel (or number of meridian)
        z0: float
            depth of the center of mass of the sphere
        clip_free_surface: bool
            if True, only mesh the part of the sphere where z < 0,
            can be used with z0 to obtain any clipped sphere.
        half: bool
            if True, only mesh the part of the sphere where y > 0
        """

        if clip_free_surface:
            if z0 < -radius: # fully immersed
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

        FloatingBody.__init__(self, nodes, panels, name=f"sphere_{next(self._ids)}")
        self.merge_duplicates()
        self.heal_triangles()


class HalfSphere(FloatingBody):
    """Floating body of the shape of half a sphere."""

    def __init__(self, **kwargs):
        Sphere.__init__(self, half=True, **kwargs)


class HorizontalCylinder(FloatingBody):
    """Floating body of the shape of a cylinder oriented along the x axis."""

    _ids = count(0)

    def __init__(self,
                 length=1.0, radius=1.0,
                 nx=11, nr=3, ntheta=11,
                 z0=0.0, clip_free_surface=False
                ):
        """Generate the mesh.

        Parameters
        ----------
        length: float
            length of the cylinder
        radius: float
            radius of the cylinder
        nx: int
            number of circular slices
        nr: int
            at the ends of the cylinder, number of points along a radius
        ntheta: int
            number of points along a circular slice of the cylinder
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
        npanels = (ntheta-1)*((nx-1)+2*max(0, (nr-1)))
        panels = np.zeros((npanels, 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta-1),
                                           range(0, nx-1))):
            panels[k, :] = (
                j+i*nx,
                j+(i+1)*nx,
                j+1+(i+1)*nx,
                j+1+i*nx
            )

        for k, (i, j) in enumerate(product(range(0, nr-1),
                                           range(ntheta*nx, ntheta*nx+ntheta-1))):
            panels[(ntheta-1)*(nx-1)+k, :] = (
                j+i*ntheta,
                j+1+i*ntheta,
                j+1+(i+1)*ntheta,
                j+(i+1)*ntheta
            )

        for k, (i, j) in enumerate(product(range(0, nr-1),
                                           range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            panels[(ntheta-1)*((nx-1)+(nr-1))+k, :] = (
                j+i*ntheta,
                j+(i+1)*ntheta,
                j+1+(i+1)*ntheta,
                j+1+i*ntheta
            )

        FloatingBody.__init__(self, nodes, panels, name=f"cylinder_{next(self._ids)}")
        self.merge_duplicates()
        self.heal_triangles()


class OneSidedRectangle(FloatingBody):
    """Rectangular panel with Cartesian mesh."""

    _ids = count(0)

    def __init__(self, height=2.0, width=10.0, nh=5, nw=5, z0=0.0):
        """Generate the mesh.

        Normals are oriented in the positive y direction.

        Parameters
        ----------
        height: float
            height of the panel (size along z)
        width: float
            width of the panel (size along x)
        nh: int
            number of points in the z direction
        nw: int
            number of points in the x direction
        z0: float
            depth of the bottom of the panel
        """

        X = np.linspace(-width/2, width/2, nw)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nw*nh, 3), dtype=np.float32)
        panels = np.zeros(((nw-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)

        FloatingBody.__init__(self, nodes, panels, name=f"rectangle_{next(self._ids)}")


class TwoSidedRectangle(FloatingBody):
    """Rectangular panel with Cartesian mesh."""

    def __init__(self, height=2.0, width=10.0, nh=5, nw=5, z0=0.0):
        """Generate the mesh.

        Parameters
        ----------
        height: float
            height of the panel (size along z)
        width: float
            width of the panel (size along x)
        nh: int
            number of points in the z direction
        nw: int
            number of points in the x direction
        z0: float
            depth of the bottom of the panel
        """

        X = np.linspace(-width/2, width/2, nw)
        Z = np.linspace(z0, z0+height, nh)

        nodes = np.zeros((nw*nh, 3), dtype=np.float32)
        panels = np.zeros((2*(nw-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)
            panels[(nw-1)*(nh-1)+k, :] = (j+i*nh, j+(i+1)*nh, j+1+(i+1)*nh, j+1+i*nh)

        FloatingBody.__init__(self, nodes, panels, name=f"rectangle_{next(OneSidedRectangle._ids)}")


class OpenRectangularParallelepiped(FloatingBody):
    """Four panels forming a parallelepiped without top nor bottom."""

    _ids = count(0)

    def __init__(self,
                 height=10.0, width=10.0, thickness=2.0,
                 nh=5, nw=5, nth=3,
                 z0=0.0
                ):
        """Generate the mesh.

        Parameters
        ----------
        height: float
            height of the object (size along z)
        width: float
            width of the object (size along x)
        thickness: float
            thickness of the object (size along y)
        nh: int
            number of points in the z direction
        nw: int
            number of points in the x direction
        nth: int
            number of points in the y direction
        z0: float
            depth of the bottom of the object
        """
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

        FloatingBody.__init__(self, combine.vertices, combine.faces, name=f"parallelepiped_{next(self._ids)}")


class RectangularParallelepiped(FloatingBody):
    """Six panels forming a complete rectangular parallelepiped."""

    def __init__(self, height=10.0, width=10.0, thickness=2.0, nh=5, nw=5, nth=3):
        """Generate the mesh.

        Parameters
        ----------
        height: float
            height of the object (size along z)
        width: float
            width of the object (size along x)
        thickness: float
            thickness of the object (size along y)
        nh: int
            number of points in the z direction
        nw: int
            number of points in the x direction
        nth: int
            number of points in the y direction
        z0: float
            depth of the bottom of the object
        """
        sides = OpenRectangularParallelepiped(height=height, width=width, thickness=thickness,
                                              nh=nh, nw=nw, nth=nth)
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

        FloatingBody.__init__(self, combine.vertices, combine.faces, name=f"parallelepiped_{next(OpenRectangularParallelepiped._ids)}")


class DummyBody(FloatingBody):
    """Body without any faces. For debugging."""

    def __init__(self):
        FloatingBody.__init__(self, np.zeros((0, 3)), np.zeros((0, 4)), name="dummy_body")
