#!/usr/bin/env python
# coding: utf-8
"""
Generate mesh for rectangles and parallelepipeds.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
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


class Rectangle(FloatingBody):
    """(One-sided) rectangle"""

    def __init__(self, size=(5.0, 5.0), resolution=(5, 5), clever=False,
                 center=(0, 0, 0), name=None):
        """Generate the mesh of a vertical rectangle.

        Normals are oriented in the positive y direction.

        Parameters
        ----------
        size : tuple of  floats
            dimensions of the rectangle (width and height)
        resolution : tuple of  ints
            number of faces along each of the two directions
        clever : bool
            if True, use the translation symmetry to speed up the computations
        center : tuple of floats
            position of the center of the rectangle
        """
        assert len(size) == 2
        assert len(center) == 3

        assert len(resolution) == 2
        assert all([isinstance(i, int) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        width, height = self.size
        self.center = np.asarray(center, dtype=np.float)
        nw, nh = resolution

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        if clever and nw > 1:
            strip = self.generate_rectangle_mesh(width=width/nw, height=height, nw=1, nh=nh, name=f"strip_of_{name}")
            strip.translate_x(-width/2 + width/(2*nw))
            mesh = TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=nw-1, name=name)
        else:
            mesh = self.generate_rectangle_mesh(width, height, nh, nw, name)
        mesh.translate(center)

        FloatingBody.__init__(self, mesh=mesh, name=name)

    def generate_rectangle_mesh(self, width=5.0, height=5.0, nh=5, nw=5, name=None):
        X = np.linspace(-width/2, width/2, nw+1)
        Z = np.linspace(-height/2, height/2, nh+1)

        nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=np.float32)
        panels = np.zeros((nw*nh, 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw), range(0, nh))):
            panels[k, :] = (j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"
        return Mesh(nodes, panels, name=f"{name}_mesh")

    @property
    def area(self):
        return self.size[0] * self.size[1]


class OpenRectangularParallelepiped(FloatingBody):
    def __init__(self, size=(5.0, 5.0, 5.0), resolution=(5, 5, 5), clever=False,
                 center=(0, 0, 0), name=None):
        """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

        Parameters
        ----------
        size : tuple of floats
            dimensions of the parallelepiped (width, thickness, height)
        resolution : tuple of ints
            number of faces along the three directions
        clever : bool
            if True, use the symmetry to speed up the computations
        center : tuple of floats
            coordinates of the center of the parallelepiped
        name : string, optional
            a name for the body
        """
        assert len(size) == 3
        assert len(center) == 3

        assert len(resolution) == 3
        assert all([isinstance(i, int) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        self.center = np.asarray(center, dtype=np.float)
        width, thickness, height = size
        nw, nth, nh = resolution

        if name is None:
            name = f"open_rectangular_parallelepiped_{next(Mesh._ids)}"

        front = Rectangle(size=(width, height), resolution=(nw, nh), clever=clever,
                          center=(0, 0, 0), name=f"front_of_{name}").mesh
        back = front.copy(name=f"back_of_{name}")

        front.translate_y(thickness/2)
        back.rotate_z(np.pi)
        back.translate_y(-thickness/2)

        if nth > 0:
            side = Rectangle(size=(thickness, height), resolution=(nth, nh), clever=clever, name=f"side_of_{name}").mesh
            other_side = side.copy(name=f"other_side_of_{name}")

            side.rotate_z(np.pi/2)
            side.translate_x(-width/2)
            other_side.rotate_z(-np.pi/2)
            other_side.translate_x(width/2)

            parallelepiped = CollectionOfMeshes([front, back, side, other_side])
        else:
            parallelepiped = CollectionOfMeshes([front, back])

        if not clever:
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)


class RectangularParallelepiped(FloatingBody):
    def __init__(self, size=(5.0, 5.0, 5.0), resolution=(5, 5, 5), clever=False,
                 center=(0, 0, 0), name=None):
        """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

        Parameters
        ----------
        size : tuple of floats
            dimensions of the parallelepiped (width, thickness, height)
        resolution : tuple of ints
            number of faces along the three directions
        clever : bool
            if True, use the symmetry to speed up the computations
        center : tuple of floats
            coordinates of the center of the parallelepiped
        name : string, optional
            a name for the body
        """
        assert len(size) == 3
        assert len(center) == 3

        assert len(resolution) == 3
        assert all([isinstance(i, int) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        self.center = np.asarray(center, dtype=np.float)
        width, thickness, height = size
        nw, nth, nh = resolution

        if name is None:
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        sides = OpenRectangularParallelepiped(size=size, resolution=resolution,
                                              clever=clever, name=f"sides_of_{name}").mesh

        if nth > 0:
            top = Rectangle(size=(width, thickness), resolution=(nw, nth),
                            clever=clever, name=f"top_of_{name}").mesh
            bottom = top.copy(name=f"bottom_of_{name}_mesh")

            top.rotate_x(np.pi/2)
            top.translate_z(height/2)
            bottom.rotate_x(-np.pi/2)
            bottom.translate_z(-height/2)

            parallelepiped = CollectionOfMeshes([sides, top, bottom])
        else:
            parallelepiped = sides

        if not clever:
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)
