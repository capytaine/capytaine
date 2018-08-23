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
from capytaine.tools.geometry import xOz_Plane, xOy_Plane, yOz_Plane

LOG = logging.getLogger(__name__)


class Rectangle(FloatingBody):
    """(One-sided) rectangle"""

    def __init__(self, size=(5.0, 5.0), resolution=(5, 5),
                 center=(0, 0, 0), normal_angles=(0, 0, 0),
                 translational_symmetry=False, reflection_symmetry=False, name=None):
        """Generate the mesh of a vertical rectangle (along y and z).

        Normals are oriented in the positive y direction.

        Parameters
        ----------
        size : couple of floats, optional
            dimensions of the rectangle (width and height)
        resolution : couple of ints, optional
            number of faces along each of the two directions
        center : 3-ple of floats, optional
            position of the center of the rectangle, default: (0, 0, 0)
        normal_angles : 3-ple of floats, optional
            direction of the normal vector, default: along x axis
        translational_symmetry : bool, optional
            if True, use the translation symmetry to speed up the computations
        reflection_symmetry : bool, optional
            if True, use the reflection symmetry to speed up the computations
        name : string, optional
            a name for the body
        """

        assert len(size) == 2, "Size of a rectangle should be given as a couple of values."
        assert all([h > 0 for h in size]), "Size of the rectangle mesh should be given as positive values."

        assert len(resolution) == 2, "Resolution of a rectangle should be given as a couple a values."
        assert all([h > 0 for h in resolution]), "Resolution of the rectangle mesh should be given as positive values."
        assert all([i == int(i) for i in resolution]), "Resolution of a rectangle should be given as integer values."

        assert len(center) == 3, "Position of the center of a rectangle should be given a 3-ple of values."

        self.size = np.asarray(size, dtype=np.float)
        width, height = self.size
        self.center = np.asarray(center, dtype=np.float)
        nw, nh = resolution

        if translational_symmetry and reflection_symmetry:
            raise NotImplemented("Rectangle generation with both reflection and translational symmetries "
                                 "has not been implemented yet.")

        if translational_symmetry and nw == 1:
            LOG.warning("To use the translation symmetry of the mesh, "
                        "it should have more than one panel in this direction. "
                        "Will return a standard mesh instead.")

        if reflection_symmetry and nw % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        if (reflection_symmetry or translational_symmetry) and normal_angles[2] != 0:
            raise ValueError("To use the symmetry of the mesh, it should be vertical.")

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        LOG.debug(f"New rectangle body of size ({width}, {height}) and resolution ({nw}, {nh}), named {name}.")

        if reflection_symmetry:
            half_mesh = Rectangle.generate_rectangle_mesh(
                width=width/2, height=height, nw=nw//2, nh=nh, center=(0, -width/4, 0), name=f"half_of_{name}_mesh"
            )
            mesh = ReflectionSymmetry(half_mesh, plane=xOz_Plane, name=f"{name}_mesh")

        elif translational_symmetry and nw > 1:
            strip = Rectangle.generate_rectangle_mesh(
                width=width/nw, height=height, nw=1, nh=nh,
                center=(0, -width/2 + width/(2*nw), 0), name=f"strip_of_{name}_mesh"
            )
            mesh = TranslationalSymmetry(strip, translation=np.asarray([0, width/nw, 0]), nb_repetitions=int(nw)-1,
                                         name=name)

        else:
            mesh = Rectangle.generate_rectangle_mesh(width=width, height=height, nw=nw, nh=nh, name=name)

        mesh.rotate(normal_angles)
        mesh.translate(center)
        FloatingBody.__init__(self, mesh=mesh, name=name)

    @staticmethod
    def generate_rectangle_mesh(width=1.0, height=1.0, nw=1, nh=1,
                                center=(0, 0, 0), normal_angles=(0, 0, 0), name=None):
        Y = np.linspace(-width/2, width/2, nw+1)
        Z = np.linspace(-height/2, height/2, nh+1)

        nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=np.float)
        panels = np.zeros((nw*nh, 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product([0.0], Y, Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw), range(0, nh))):
            panels[k, :] = (j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        mesh = Mesh(nodes, panels, name=f"{name}_mesh")
        mesh.rotate(normal_angles)
        mesh.translate(center)

        return mesh

    @property
    def area(self):
        return self.size[0] * self.size[1]


class RectangularParallelepiped(FloatingBody):
    def __init__(self,
                 size=(5.0, 5.0, 5.0), resolution=(5, 5, 5),
                 center=(0, 0, 0),
                 top=True, bottom=True,
                 reflection_symmetry=False,
                 translational_symmetry=False,
                 name=None):
        """Generate the mesh of six rectangles forming a parallelepiped.

        Parameters
        ----------
        size : 3-ple of floats, optional
            dimensions of the parallelepiped (width, thickness, height) for coordinates (x, y, z).
        resolution : 3-ple of ints, optional
            number of faces along the three directions
        center : 3-ple of floats, optional
            coordinates of the center of the parallelepiped
        top: bool, optional
            whether or not to close the parallelepiped on the top
        bottom: bool, optional
            whether or not to close the parallelepiped on the bottom
        reflection_symmetry : bool, optional
            use xOz and yOz symmetry plane to generate the mesh
        translational_symmetry : bool, optional
            if True, use the translation symmetry in the x direction to speed up the computations.
            To use the translation symmetry in the y direction, create a x-symmetric body and then rotate it by pi/2.
        name : string, optional
            a name for the body
        """

        assert len(size) == 3, "Size of a rectangular parallelepiped should be given as a 3-ple of values."
        assert all([h > 0 for h in size]), "Size of the rectangular mesh should be given as positive values."

        assert len(resolution) == 3, "Resolution of a rectangular parallelepiped should be given as a 3-ple a values."
        assert all([h > 0 for h in resolution]), "Resolution of the rectangular parallelepiped mesh " \
                                                 "should be given as positive values."
        assert all([i == int(i) for i in resolution]), "Resolution of a rectangular parallelepiped " \
                                                       "should be given as integer values."

        assert len(center) == 3, "Position of the center of a parallelepiped should be given a 3-ple of values."

        self.size = np.asarray(size, dtype=np.float)
        width, thickness, height = size
        self.center = np.asarray(center, dtype=np.float)
        nw, nth, nh = resolution

        if translational_symmetry and reflection_symmetry:
            raise NotImplemented("Parallelepiped generation with both reflection and translational symmetries "
                                 "has not been implemented yet.")

        if reflection_symmetry:
            raise NotImplemented("TODO")

        if reflection_symmetry and (nw % 2 == 1 or nth % 2 == 1):
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        if name is None:
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        LOG.debug(f"New rectangular parallelepiped body "
                  f"of size ({width}, {thickness}, {height}) and resolution ({resolution}), named {name}.")

        front_panel = Rectangle.generate_rectangle_mesh(
            width=width/nw, height=height, nw=1, nh=nh,
            center=(-width/2 + width/(2*nw), thickness/2, 0),
            normal_angles=(0, 0, -np.pi/2),
            name=f"front_panel_of_{name}_mesh"
        )

        back_panel = front_panel.mirror(plane=xOz_Plane, inplace=False, name=f"back_panel_of_{name}_mesh")

        top_panel = Rectangle.generate_rectangle_mesh(
            width=thickness, height=width/nw, nw=nth, nh=1,
            center=(-width/2 + width/(2*nw), 0, height/2),
            normal_angles=(0, np.pi/2, 0),
            name=f"top_panel_of_{name}_mesh"
        )

        bottom_panel = top_panel.mirror(plane=xOy_Plane, inplace=False, name=f"bottom_panel_of_{name}_mesh")

        panels = [front_panel, back_panel]
        if top:
            panels.append(top_panel)
        if bottom:
            panels.append(bottom_panel)
        four_panels = CollectionOfMeshes(panels, name=f"ring_of_{name}_mesh")

        front_back_top_bottom = TranslationalSymmetry(
            four_panels,
            translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1,
            name=f"body_of_{name}_mesh"
        )

        side = Rectangle.generate_rectangle_mesh(
            width=thickness, height=height, nw=nth, nh=nh,
            center=(width/2, 0, 0),
            name=f"side_of_{name}_mesh"
        )
        other_side = side.mirror(plane=yOz_Plane, inplace=False, name=f"other_side_of_{name}_mesh")

        parallelepiped = CollectionOfMeshes([front_back_top_bottom, side, other_side], name=f"{name}_mesh")

        if not (reflection_symmetry or translational_symmetry):
            parallelepiped = parallelepiped.merge(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    @property
    def volume(self):
        return np.product(self.size)


class OpenRectangularParallelepiped(RectangularParallelepiped):
    def __init__(self, *args, **kwargs):
        RectangularParallelepiped.__init__(self, top=False, bottom=False, *args, **kwargs)
        # Kept mostly for legacy

