#!/usr/bin/env python
# coding: utf-8
"""Generate mesh for rectangles and parallelepipeds."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging
from itertools import product

import numpy as np

from capytaine.mesh.mesh import Mesh
from capytaine.mesh.meshes_collection import CollectionOfMeshes
from capytaine.mesh.symmetries import TranslationalSymmetry, ReflectionSymmetry
from capytaine.bodies import FloatingBody
from capytaine.tools.geometry import xOz_Plane, yOz_Plane

LOG = logging.getLogger(__name__)


class Rectangle(FloatingBody):
    """(One-sided) rectangle"""

    def __init__(self, size=(5.0, 5.0), resolution=(5, 5), center=(0, 0, 0),
                 clever=False, name=None):
        """Generate the mesh of a vertical rectangle (along x and z).

        Normals are oriented in the positive y direction.

        Parameters
        ----------
        size : tuple of  floats, optional
            dimensions of the rectangle (width and height)
        resolution : tuple of  ints, optional
            number of faces along each of the two directions
        center : tuple of floats, optional
            position of the center of the rectangle
        clever : bool, optional
            if True, use the translation symmetry along the x axis to speed up the computations
        name : string, optional
            a name for the body
        """
        assert len(size) == 2
        assert len(center) == 3

        assert len(resolution) == 2
        assert all([i == int(i) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        width, height = self.size
        self.center = np.asarray(center, dtype=np.float)
        nw, nh = resolution

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        if clever and nw > 1:
            strip = self.generate_rectangle_mesh(width=width/nw, height=height, nw=1, nh=nh, name=f"strip_of_{name}")
            strip.translate_x(-width/2 + width/(2*nw))
            mesh = TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=int(nw)-1, name=name)
        else:
            mesh = self.generate_rectangle_mesh(width, height, nh, nw, name)
        mesh.translate(center)

        FloatingBody.__init__(self, mesh=mesh, name=name)

    def generate_rectangle_mesh(self, width=5.0, height=5.0, nh=5, nw=5, name=None):
        X = np.linspace(-width/2, width/2, nw+1)
        Z = np.linspace(-height/2, height/2, nh+1)

        nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=np.float)
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
    def __init__(self, size=(5.0, 5.0, 5.0), resolution=(5, 5, 5), center=(0, 0, 0),
                 clever=False, name=None):
        """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

        Parameters
        ----------
        size : tuple of floats, optional
            dimensions of the parallelepiped (width, thickness, height) or (dx, dy, dz)
        resolution : tuple of ints, optional
            number of faces along the three directions
        center : tuple of floats, optional
            coordinates of the center of the parallelepiped
        clever : bool, optional
            if True, use the translation symmetry in the x direction to speed up the computations
            To use the translation symmetry in the y direction, create a x-symmetric body and then rotate it by pi/2.
        name : string, optional
            a name for the body
        """
        assert len(size) == 3
        assert len(center) == 3

        assert len(resolution) == 3
        assert all([i == int(i) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        self.center = np.asarray(center, dtype=np.float)
        width, thickness, height = size
        nw, nth, nh = resolution

        if name is None:
            name = f"open_rectangular_parallelepiped_{next(Mesh._ids)}"

        front_panel = Rectangle(size=(width/nw, height), resolution=(1, nh), center=(0, 0, 0),
                                clever=False, name=f"front_panel_of_{name}").mesh
        back_panel = front_panel.copy(name=f"back_panel_of_{name}_mesh")
        front_panel.translate((-width/2 + width/(2*nw), thickness/2, 0))
        back_panel.rotate_z(np.pi)
        back_panel.translate((-width/2 + width/(2*nw), -thickness/2, 0))
        front_and_back_panel = CollectionOfMeshes([front_panel, back_panel], name="front_and_back_panels_of_{name}_mesh")

        front_and_back = TranslationalSymmetry(front_and_back_panel,
                                               translation=(width/nw, 0, 0),
                                               nb_repetitions=int(nw)-1,
                                               name=f"front_and_back_of_{name}_mesh")

        if nth > 0:
            side = Rectangle(size=(thickness, height), resolution=(nth, nh), center=(0, 0, 0),
                             clever=False, name=f"side_of_{name}").mesh
            other_side = side.copy(name=f"other_side_of_{name}_mesh")

            side.rotate_z(np.pi/2)
            side.translate_x(-width/2)
            other_side.rotate_z(-np.pi/2)
            other_side.translate_x(width/2)

            parallelepiped = CollectionOfMeshes([front_and_back, side, other_side])
        else:
            parallelepiped = CollectionOfMeshes([front_and_back])

        if not clever:
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        parallelepiped.name = f"{name}_mesh"
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)


class RectangularParallelepiped(FloatingBody):
    def __init__(self, size=(5.0, 5.0, 5.0), resolution=(5, 5, 5), center=(0, 0, 0),
                 clever=False, name=None):
        """Generate the mesh of four panels forming a parallelepiped without top nor bottom.

        Parameters
        ----------
        size : tuple of floats
            dimensions of the parallelepiped (width, thickness, height)
        resolution : tuple of ints
            number of faces along the three directions
        center : tuple of floats
            coordinates of the center of the parallelepiped
        clever : bool
            if True, use the translation symmetry in the x direction to speed up the computations
            To use the translation symmetry in the y direction, create a x-symmetric body and then rotate it by pi/2.
        name : string, optional
            a name for the body
        """
        assert len(size) == 3
        assert len(center) == 3

        assert len(resolution) == 3
        assert all([i == int(i) for i in resolution])

        self.size = np.asarray(size, dtype=np.float)
        self.center = np.asarray(center, dtype=np.float)
        width, thickness, height = size
        nw, nth, nh = resolution

        if name is None:
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        front_panel = Rectangle(size=(width/nw, height), resolution=(1, nh), center=(0, 0, 0),
                                clever=False, name=f"front_panel_of_{name}").mesh
        back_panel = front_panel.copy(name=f"back_panel_of_{name}_mesh")
        front_panel.translate((-width/2 + width/(2*nw), thickness/2, 0))
        back_panel.rotate_z(np.pi)
        back_panel.translate((-width/2 + width/(2*nw), -thickness/2, 0))

        if nth > 0:
            top = Rectangle(size=(width/nw, thickness), resolution=(1, nth),
                            clever=False, name=f"top_panel_of_{name}").mesh
            bottom = top.copy(name=f"bottom_panel_of_{name}_mesh")

            top.translate_x(-width/2 + width/(2*nw))
            top.rotate_x(np.pi/2)
            top.translate_z(height/2)
            bottom.translate_x(-width/2 + width/(2*nw))
            bottom.rotate_x(-np.pi/2)
            bottom.translate_z(-height/2)

            four_panels = CollectionOfMeshes([front_panel, back_panel, top, bottom], name=f"ring_of_{name}_mesh")

            front_back_top_bottom = TranslationalSymmetry(four_panels,
                                                          translation=(width/nw, 0, 0),
                                                          nb_repetitions=int(nw)-1,
                                                          name=f"body_of_{name}_mesh")

            side = Rectangle(size=(thickness, height), resolution=(nth, nh), center=(0, 0, 0),
                             clever=False, name=f"side_of_{name}").mesh
            other_side = side.copy(name=f"other_side_of_{name}_mesh")

            side.rotate_z(np.pi/2)
            side.translate_x(-width/2)
            other_side.rotate_z(-np.pi/2)
            other_side.translate_x(width/2)

            parallelepiped = CollectionOfMeshes([front_back_top_bottom, side, other_side])
        else:
            front_and_back = TranslationalSymmetry(front_panel + back_panel, translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1)
            parallelepiped = CollectionOfMeshes([front_and_back])

        if not clever:
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        parallelepiped.name = f"{name}_mesh"
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)
