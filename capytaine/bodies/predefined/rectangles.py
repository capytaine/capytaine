"""Legacy interfaces to predefined meshes"""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
import numpy as np

from capytaine.meshes.predefined import mesh_rectangle, mesh_parallelepiped
from capytaine.meshes.meshes import Mesh
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

        LOG.warning("Deprecation warning: The class Rectangle() is deprecated. "
                    "Please prefer the function capytaine.meshes.predefined.mesh_rectangle()")

        self.size = np.asarray(size, dtype=float)
        self.geometric_center = np.asarray(center, dtype=float)

        if name is None:
            name = f"rectangle_{next(Mesh._ids)}"

        mesh = mesh_rectangle(size=size, resolution=resolution, center=center, normal=normal,
                translation_symmetry=translational_symmetry, reflection_symmetry=reflection_symmetry,
                name=f"{name}_mesh")
        FloatingBody.__init__(self, mesh=mesh, name=name)


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

        LOG.warning("Deprecation warning: The class RectangularParallelepiped() is deprecated. "
                    "Please prefer the function capytaine.meshes.predefined.mesh_parallelepiped()")

        if name is None:
            name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

        missing_sides = set()
        if not top: missing_sides.add("top")
        if not bottom: missing_sides.add("bottom")

        mesh = mesh_parallelepiped(size=size, resolution=resolution, center=center,
                missing_sides=missing_sides,
                translation_symmetry=translational_symmetry, reflection_symmetry=reflection_symmetry,
                name=f"{name}_mesh")

        FloatingBody.__init__(self, mesh=mesh, name=name)

class OpenRectangularParallelepiped(RectangularParallelepiped):
    def __init__(self, *args, **kwargs):
        RectangularParallelepiped.__init__(self, top=False, bottom=False, *args, **kwargs)
        # Kept mostly for legacy
