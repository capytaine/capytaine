"""Legacy interfaces to predefined meshes"""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
import numpy as np

from capytaine.meshes.predefined import mesh_disk, mesh_vertical_cylinder, mesh_horizontal_cylinder
from capytaine.meshes.meshes import Mesh
from capytaine.bodies.bodies import FloatingBody

LOG = logging.getLogger(__name__)


##########
#  Disk  #
##########

class Disk(FloatingBody):
    """(One-sided) disk.
    Deprecated: please prefer capytaine.meshes.predefined.mesh_disk()

    Parameters
    ----------
    radius : float, optional
        radius of the disk
    resolution : 2-ple of int, optional
        number of panels along a radius and around the disk
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the disk
    normal: 3-ple of floats, optional
        normal vector, default: along x axis
    axial_symmetry : bool, optional
        if True, use the axial symmetry to speed up the computations
    reflection_symmetry : bool, optional
        if True, use the reflection symmetry to speed up the computations
    name : str, optional
        a string naming the floating body
    """

    def __init__(self, radius=1.0, resolution=(3, 5),
                 center=(0, 0, 0), normal=(1, 0, 0),
                 reflection_symmetry=False, axial_symmetry=False,
                 name=None):
        LOG.warning("Deprecation warning: The class Disk() is deprecated. "
                "Please prefer the function capytaine.meshes.predefined.mesh_disk()")

        if name is None:
            name = f"disk_{next(Mesh._ids)}"

        self.radius = float(radius)
        self.geometric_center = np.asarray(center, dtype=float)
        mesh = mesh_disk(radius=radius, center=center, normal=normal, resolution=resolution,
                reflection_symmetry=reflection_symmetry, axial_symmetry=axial_symmetry, name=f"{name}_mesh")
        FloatingBody.__init__(self, mesh=mesh, name=name)


##############
#  Cylinder  #
##############

class HorizontalCylinder(FloatingBody):
    """Horizontal cylinder
    Deprecated: please prefer capytaine.meshes.predefined.mesh_horizontal_cylinder()

    Parameters
    ----------
    length : float, optional
        length of the cylinder
    radius : float, optional
        radius of the cylinder
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the cylinder
    nx : int, optional
        number of circular slices
    ntheta : int, optional
        number of panels along a circular slice of the cylinder
    nr : int, optional
        number of panels along a radius on the extremities of the cylinder
    reflection_symmetry : bool, optional
        if True, returns a ReflectionSymmetricMesh
    translation_symmetry : bool, optional
        if True, uses a TranslationalSymmetricMesh internally for the main part of the cylinder
    name : str, optional
        a string naming the floating body
    """

    def __init__(self, length=10.0, radius=1.0, center=(0, 0, 0),
                 nx=10, ntheta=10, nr=2,
                 reflection_symmetry=True, translation_symmetry=False,
                 clever=None,
                 name=None):

        LOG.warning("Deprecation warning: The class HorizontalCylinder() is deprecated. "
                "Please prefer the function capytaine.meshes.predefined.mesh_horizontal_cylinder()")

        self.length = length
        self.radius = radius
        self.geometric_center = np.asarray(center, dtype=float)

        if name is None:
            name = f"cylinder_{next(Mesh._ids)}"

        mesh = mesh_horizontal_cylinder(length=length, radius=radius, center=center,
                resolution=(nr, ntheta, nx), reflection_symmetry=reflection_symmetry,
                translation_symmetry=translation_symmetry, name=f"{name}_mesh")
        FloatingBody.__init__(self, mesh=mesh, name=name)


class VerticalCylinder(FloatingBody):
    """Vertical cylinder.
    Deprecated: please prefer capytaine.meshes.predefined.mesh_vertical_cylinder()

    Parameters
    ----------
    length : float, optional
        length of the cylinder
    radius : float, optional
        radius of the cylinder
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the cylinder
    nx : int, optional
        number of circular slices
    ntheta : int, optional
        number of panels along a circular slice of the cylinder
    nr : int, optional
        number of panels along a radius on the extremities of the cylinder
    clever : bool, optional
        if True, uses the mesh symmetries
    name : str, optional
        a string naming the floating body
    """

    def __init__(self, length=10.0, radius=1.0, center=(0, 0, 0),
                 nx=10, ntheta=10, nr=2,
                 clever=True, name=None):
        LOG.warning("Deprecation warning: The class VerticalCylinder() is deprecated. "
                "Please prefer the function capytaine.meshes.predefined.mesh_vertical_cylinder()")

        self.length = length
        self.radius = radius
        self.geometric_center = np.asarray(center, dtype=float)

        if name is None:
            name = f"cylinder_{next(Mesh._ids)}"

        mesh = mesh_vertical_cylinder(length=length, radius=radius, center=center,
                resolution=(nr, ntheta, nx), reflection_symmetry=False,
                axial_symmetry=clever, name=f"{name}_mesh")

        FloatingBody.__init__(self, mesh=mesh, name=name)
