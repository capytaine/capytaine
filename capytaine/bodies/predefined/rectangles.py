#!/usr/bin/env python
# coding: utf-8
"""Generate mesh of rectangles and parallelepipeds."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from itertools import product

import numpy as np

from capytaine.meshes.geometry import xOz_Plane, xOy_Plane, yOz_Plane, e_x, Ox_axis
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import TranslationalSymmetricMesh, ReflectionSymmetricMesh
from capytaine.bodies.bodies import FloatingBody

LOG = logging.getLogger(__name__)


class Rectangle(FloatingBody):
    """One-sided vertical rectangle (along y and z).

    By default, the normals are oriented in the positive y direction.

    Parameters
    ----------
    size : couple of floats, optional
        dimensions of the rectangle (width and height)
    resolution : couple of ints, optional
        number of faces along each of the two directions
    center : 3-ple of floats, optional
        position of the geometric center of the rectangle, default: (0, 0, 0)
    normal: 3-ple of floats, optional
        normal vector, default: along x axis
    translational_symmetry : bool, optional
        if True, use the translation symmetry to speed up the computations
    reflection_symmetry : bool, optional
        if True, use the reflection symmetry to speed up the computations
    name : string, optional
        a name for the body
    """

    def __init__(self, size=(5.0, 5.0), resolution=(5, 5),
                 center=(0, 0, 0), normal=(1, 0, 0),
                 translational_symmetry=False, reflection_symmetry=False, name=None):
        assert len(size) == 2, "Size of a rectangle should be given as a couple of values."
        assert all([h > 0 for h in size]), "Size of the rectangle mesh should be given as positive values."

        assert len(resolution) == 2, "Resolution of a rectangle should be given as a couple a values."
        assert all([h > 0 for h in resolution]), "Resolution of the rectangle mesh should be given as positive values."
        assert all([i == int(i) for i in resolution]), "Resolution of a rectangle should be given as integer values."

        assert len(center) == 3, "Position of the center of a rectangle should be given a 3-ple of values."

        self.size = np.asarray(size, dtype=float)
        width, height = self.size
        self.geometric_center = np.asarray(center, dtype=float)
        nw, nh = resolution

        if translational_symmetry and reflection_symmetry:
            raise NotImplementedError("Rectangle generation with both reflection and translational symmetries "
                                      "has not been implemented yet.")

        if translational_symmetry and nw == 1:
            LOG.warning("To use the translation symmetry of the mesh, "
                        "it should have more than one panel in this direction. "
                        "Will return a standard mesh instead.")

        if reflection_symmetry and nw % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        if (reflection_symmetry or translational_symmetry) and normal[2] != 0:
            raise ValueError("To use the symmetry of the mesh, it should be vertical.")

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        LOG.debug(f"New rectangle body of size ({width}, {height}) and resolution ({nw}, {nh}), named {name}.")

        if reflection_symmetry:
            half_mesh = Rectangle.generate_rectangle_mesh(
                width=width/2, height=height, nw=nw//2, nh=nh,
                center=(0, -width/4, 0), name=f"half_of_{name}_mesh"
            )
            mesh = ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=f"{name}_mesh")

        elif translational_symmetry and nw > 1:
            strip = Rectangle.generate_rectangle_mesh(
                width=width/nw, height=height, nw=1, nh=nh,
                center=(0, -width/2 + width/(2*nw), 0), name=f"strip_of_{name}_mesh"
            )
            mesh = TranslationalSymmetricMesh(strip,
                                              translation=np.asarray([0, width/nw, 0]), nb_repetitions=int(nw)-1,
                                              name=name)

        else:
            mesh = Rectangle.generate_rectangle_mesh(width=width, height=height, nw=nw, nh=nh, name=name)

        mesh.rotate_around_center_to_align_vectors((0, 0, 0), mesh.faces_normals[0], normal)
        mesh.translate(center)
        FloatingBody.__init__(self, mesh=mesh, name=name)

    @staticmethod
    def generate_rectangle_mesh(width=1.0, height=1.0, nw=1, nh=1,
                                center=(0, 0, 0), normal=(1, 0, 0), name=None):
        Y = np.linspace(-width/2, width/2, nw+1)
        Z = np.linspace(-height/2, height/2, nh+1)

        nodes = np.zeros(((nw+1)*(nh+1), 3), dtype=float)
        panels = np.zeros((nw*nh, 4), dtype=int)

        for i, (x, y, z) in enumerate(product([0.0], Y, Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nw), range(0, nh))):
            panels[k, :] = (j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        mesh = Mesh(nodes, panels, name=f"{name}_mesh")
        mesh.rotate_around_center_to_align_vectors((0, 0, 0), mesh.faces_normals[0], normal)
        mesh.translate(center)

        return mesh


class RectangularParallelepiped(FloatingBody):
    """Six rectangles forming a parallelepiped.

    Parameters
    ----------
    size : 3-ple of floats, optional
        dimensions of the parallelepiped (width, thickness, height) for coordinates (x, y, z).
    resolution : 3-ple of ints, optional
        number of faces along the three directions
    center : 3-ple of floats, optional
        coordinates of the geometric center of the parallelepiped
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

    def __init__(self,
                 size=(1.0, 1.0, 1.0), resolution=(4, 4, 4),
                 center=(0, 0, 0),
                 top=True, bottom=True,
                 reflection_symmetry=False,
                 translational_symmetry=False,
                 name=None):

        assert len(size) == 3, "Size of a rectangular parallelepiped should be given as a 3-ple of values."
        assert all([h > 0 for h in size]), "Size of the rectangular mesh should be given as positive values."

        assert len(resolution) == 3, "Resolution of a rectangular parallelepiped should be given as a 3-ple a values."
        assert all([h > 0 for h in resolution]), "Resolution of the rectangular parallelepiped mesh " \
                                                 "should be given as positive values."
        assert all([i == int(i) for i in resolution]), "Resolution of a rectangular parallelepiped " \
                                                       "should be given as integer values."

        assert len(center) == 3, "Position of the center of a parallelepiped should be given a 3-ple of values."

        self.size = np.asarray(size, dtype=float)
        width, thickness, height = size
        self.geometric_center = np.asarray(center, dtype=float)
        nw, nth, nh = resolution

        if translational_symmetry and reflection_symmetry:
            raise NotImplementedError("Parallelepiped generation with both reflection and translational symmetries "
                                      "has not been implemented yet.")

        if reflection_symmetry and (nw % 2 == 1 or nth % 2 == 1):
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        if name is None:
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        LOG.debug(f"New rectangular parallelepiped body "
                  f"of size ({width}, {thickness}, {height}) and resolution ({resolution}), named {name}.")

        if reflection_symmetry:
            parallelepiped = self._generate_mesh_with_reflection_symmetry(resolution, top, bottom, name)
        else:
            parallelepiped = self._generate_mesh_with_translational_symmetry(resolution, top, bottom, name)

        if not (reflection_symmetry or translational_symmetry):
            parallelepiped = parallelepiped.merged(name=f"{name}_mesh")
            parallelepiped.merge_duplicates()
            parallelepiped.heal_triangles()

        parallelepiped.translate(center)
        FloatingBody.__init__(self, mesh=parallelepiped, name=name)

    def _generate_mesh_with_translational_symmetry(self, resolution, top, bottom, name):
        width, thickness, height = self.size
        nw, nth, nh = resolution

        front_panel = Rectangle.generate_rectangle_mesh(
            width=width/nw, height=height, nw=1, nh=nh,
            center=(-width/2 + width/(2*nw), thickness/2, 0),
            normal=(0, 1, 0),
            name=f"front_panel_of_{name}_mesh"
        )

        back_panel = front_panel.mirror(plane=xOz_Plane, inplace=False, name=f"back_panel_of_{name}_mesh")

        top_panel = Rectangle.generate_rectangle_mesh(
            width=thickness, height=width/nw, nw=nth, nh=1,
            center=(-width/2 + width/(2*nw), 0, height/2),
            normal=(0, 0, 1),
            name=f"top_panel_of_{name}_mesh"
        )

        bottom_panel = top_panel.mirror(plane=xOy_Plane, inplace=False, name=f"bottom_panel_of_{name}_mesh")

        panels = [front_panel, back_panel]
        if top:
            panels.append(top_panel)
        if bottom:
            panels.append(bottom_panel)
        ring = CollectionOfMeshes(panels, name=f"ring_of_{name}_mesh").merged()
        ring.merge_duplicates()
        ring.heal_triangles()

        open_parallelepiped = TranslationalSymmetricMesh(
            ring,
            translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1,
            name=f"body_of_{name}_mesh"
        )

        side = Rectangle.generate_rectangle_mesh(
            width=thickness, height=height, nw=nth, nh=nh,
            center=(width/2, 0, 0),
            normal=(1, 0, 0),
            name=f"side_of_{name}_mesh"
        )

        other_side = side.mirror(plane=yOz_Plane, inplace=False, name=f"other_side_of_{name}_mesh")

        return CollectionOfMeshes([open_parallelepiped, side, other_side], name=f"{name}_mesh")

    def _generate_mesh_with_reflection_symmetry(self, resolution, top, bottom, name):
        width, thickness, height = self.size
        nw, nth, nh = resolution

        half_front = Rectangle.generate_rectangle_mesh(
            width=width/2, height=height, nw=nw//2, nh=nh,
            center=(-width/4, thickness/2, 0),
            normal=(0, 1, 0),
            name=f"half_front_of_{name}_mesh"
        )

        quarter_of_top = Rectangle.generate_rectangle_mesh(
            width=thickness/2, height=width/2, nw=nth//2, nh=nw//2,
            center=(-width/4, thickness/4, height/2),
            normal=(0, 0, 1),
            name=f"top_panel_of_{name}_mesh"
        )

        quarter_of_bottom = quarter_of_top.mirror(plane=xOy_Plane, inplace=False, name=f"bottom_panel_of_{name}_mesh")

        half_side = Rectangle.generate_rectangle_mesh(
            width=thickness/2, height=height, nw=nth//2, nh=nh,
            center=(-width/2, thickness/4, 0),
            normal=(1, 0, 0),
            name=f"half_side_of_{name}_mesh"
        )

        panels = [half_front, half_side]
        if top:
            panels.append(quarter_of_top)
        if bottom:
            panels.append(quarter_of_bottom)
        quarter_of_mesh = CollectionOfMeshes(panels, name=f"quarter_of_{name}_mesh").merged()

        half_mesh = ReflectionSymmetricMesh(quarter_of_mesh, plane=yOz_Plane, name=f"half_of_{name}_mesh")
        return ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=f"{name}_mesh")


class OpenRectangularParallelepiped(RectangularParallelepiped):
    def __init__(self, *args, **kwargs):
        RectangularParallelepiped.__init__(self, top=False, bottom=False, *args, **kwargs)
        # Kept mostly for legacy
