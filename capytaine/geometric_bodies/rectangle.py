#!/usr/bin/env python
# coding: utf-8
"""Generate mesh of rectangles and parallelepipeds."""
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
                 translation_symmetry=False, reflection_symmetry=False, name=None):
        """Generate the mesh of a vertical rectangle (along x and z).

        Normals are oriented in the positive y direction.

        Parameters
        ----------
        size : couple of floats, optional
            dimensions of the rectangle (width and height)
        resolution : couple of ints, optional
            number of faces along each of the two directions
        center : 3-ple of floats, optional
            position of the center of the rectangle
        translation_symmetry : bool, optional
            if True, use the translation symmetry along the x axis to speed up the computations
        reflection_symmetry : bool, optional
            if True, use the reflection symmetry across the yOz plane to speed up the computations
        name : string, optional
            a name for the body
        """

        assert len(size) == 2, "Size of a rectangle should be given as a couple of values."
        assert all([h > 0 for h in size]), "Size of the rectangular mesh should be given as positive values."

        assert len(resolution) == 2, "Resolution of a rectangle should be given as a couple a values."
        assert all([h > 0 for h in size]), "Resolution of the rectangular mesh should be given as positive values."
        assert all([i == int(i) for i in resolution]), "Resolution of a rectangle should be given as integer values."

        assert len(center) == 3, "Position of the center of a rectangle should be given a 3-ple of values."

        self.size = np.asarray(size, dtype=np.float)
        width, height = self.size
        self.center = np.asarray(center, dtype=np.float)
        nw, nh = resolution

        if translation_symmetry and reflection_symmetry:
            raise NotImplemented("Rectangle generation with both reflection and translation symmetries "
                                 "has not been implemented yet.")

        if translation_symmetry and nw == 1:
            LOG.warning("To use the translation symmetry of the mesh, "
                        "it should have more than one panel in this direction. "
                        "Will return a standard mesh instead.")

        if reflection_symmetry and nw % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        LOG.debug(f"New rectangular body of size ({width}, {height}) and resolution ({nw}, {nh}), named {name}.")

        if reflection_symmetry:
            half_mesh = Rectangle.generate_rectangle_mesh(
                width=width/2, height=height, nw=nw//2, nh=nh, name=f"half_of_{name}_mesh"
            )
            half_mesh.translate_x(-width/4)
            mesh = ReflectionSymmetry(half_mesh, plane=yOz_Plane, name=f"{name}_mesh")

        elif translation_symmetry and nw > 1:
            strip = Rectangle.generate_rectangle_mesh(
                width=width/nw, height=height, nw=1, nh=nh, name=f"strip_of_{name}_mesh"
            )
            strip.translate_x(-width/2 + width/(2*nw))
            mesh = TranslationalSymmetry(strip, translation=np.asarray([width/nw, 0.0, 0.0]), nb_repetitions=int(nw)-1,
                                         name=name)

        else:
            mesh = Rectangle.generate_rectangle_mesh(width=width, height=height, nw=nw, nh=nh, name=name)

        mesh.translate(center)
        FloatingBody.__init__(self, mesh=mesh, name=name)

    @staticmethod
    def generate_rectangle_mesh(width=5.0, height=5.0, nw=5, nh=5, name=None):
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
    def __init__(self, *args, **kwargs):
        RectangularParallelepiped.__init__(self, top=False, bottom=False, *args, **kwargs)
        # Kept mostly for legacy


class RectangularParallelepiped(FloatingBody):
    def __init__(self,
                 size=(5.0, 5.0, 5.0), resolution=(5, 5, 5),
                 center=(0, 0, 0),
                 top=True, bottom=True,
                 clever=False, name=None):
        """Generate the mesh of six panels forming a parallelepiped

        Parameters
        ----------
        size : tuple of floats, optional
            dimensions of the parallelepiped (width, thickness, height)
        resolution : tuple of ints, optional
            number of faces along the three directions
        center : tuple of floats, optional
            coordinates of the center of the parallelepiped
        top: bool, optional
            whether or not to close the parallelepiped on the top
        bottom: bool, optional
            whether or not to close the parallelepiped on the bottom
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
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        front_panel = Rectangle.generate_rectangle_mesh(
            width=width/nw, height=height, nw=1, nh=nh, name=f"front_panel_of_{name}_mesh"
        )
        front_panel.translate((-width/2 + width/(2*nw), thickness/2, 0))

        back_panel = Rectangle.generate_rectangle_mesh(
            width=width/nw, height=height, nw=1, nh=nh, name=f"back_panel_of_{name}_mesh"
        )
        back_panel.rotate_z(np.pi)
        back_panel.translate((-width/2 + width/(2*nw), -thickness/2, 0))

        panels = [front_panel, back_panel]

        if nth > 0:

            if top:
                top_panel = Rectangle.generate_rectangle_mesh(
                    width=width/nw, height=thickness, nw=1, nh=nth, name=f"top_panel_of_{name}_mesh"
                )
                top_panel.translate_x(-width/2 + width/(2*nw))
                top_panel.rotate_x(np.pi/2)
                top_panel.translate_z(height/2)
                panels.append(top_panel)

            if bottom:
                bottom_panel = Rectangle.generate_rectangle_mesh(
                    width=width/nw, height=thickness, nw=1, nh=nth, name=f"bottom_panel_of_{name}_mesh"
                )
                bottom_panel.translate_x(-width/2 + width/(2*nw))
                bottom_panel.rotate_x(-np.pi/2)
                bottom_panel.translate_z(-height/2)
                panels.append(bottom_panel)

            four_panels = CollectionOfMeshes(panels, name=f"ring_of_{name}_mesh")

            front_back_top_bottom = TranslationalSymmetry(
                four_panels,
                translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1,
                name=f"body_of_{name}_mesh"
            )

            side = Rectangle.generate_rectangle_mesh(
                width=thickness, height=height, nw=nth, nh=nh, name=f"side_of_{name}_mesh"
            )
            other_side = side.copy(name=f"other_side_of_{name}_mesh")

            side.rotate_z(np.pi/2)
            side.translate_x(-width/2)
            other_side.rotate_z(-np.pi/2)
            other_side.translate_x(width/2)

            parallelepiped = CollectionOfMeshes([front_back_top_bottom, side, other_side], name=f"{name}_mesh")
        else:
            front_and_back = TranslationalSymmetry(front_panel + back_panel, translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1)
            parallelepiped = CollectionOfMeshes([front_and_back], name=f"{name}_mesh")

        if not clever:
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)
