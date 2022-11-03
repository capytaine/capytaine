"""Generate meshes of cylinders and disks"""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
from itertools import product

import numpy as np
from numpy import pi, cos, sin

from capytaine.meshes.geometry import xOy_Plane, xOz_Plane, yOz_Plane, Oz_axis
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import TranslationalSymmetricMesh, AxialSymmetricMesh, ReflectionSymmetricMesh

LOG = logging.getLogger(__name__)


def mesh_disk(*, radius=1.0, center=(0, 0, 0), normal=(0, 0, 1), resolution=(3, 6),
        reflection_symmetry=False, axial_symmetry=False, name=None, _theta_max=2*pi):
    """(One-sided) disk.

    Parameters
    ----------
    radius : float, optional
        radius of the disk
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the disk
    normal: 3-ple of floats, optional
        normal vector, default: along x axis
    resolution : 2-ple of int, optional
        number of panels along a radius and around the disk
    axial_symmetry : bool, optional
        if True, returns an AxialSymmetricMesh
    reflection_symmetry : bool, optional
        if True, returns a ReflectionSymmetricMesh
    name : str, optional
        a string naming the mesh
    _theta_max: float, optional
        internal parameter, to return an arc circle instead of a full circle
    """
    assert radius > 0, "Radius of the disk mesh should be given as a positive value."

    assert len(resolution) == 2, "Resolution of a disk should be given as a couple of values."
    assert all([h > 0 for h in resolution]), "Resolution of the disk mesh should be given as positive values."
    assert all([i == int(i) for i in resolution]), "Resolution of a disk should be given as integer values."

    assert len(center) == 3, "Position of the center of a disk should be given a 3-ple of values."

    nr, ntheta = resolution

    if name is None:
        name = f"disk_{next(Mesh._ids)}"

    if reflection_symmetry and axial_symmetry:
        raise NotImplementedError("Disks with both symmetries have not been implemented.")

    LOG.debug(f"New disk of radius {radius} and resolution {resolution}, named {name}.")

    if reflection_symmetry:
        if ntheta % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        half_mesh = mesh_disk(radius=radius, _theta_max=_theta_max/2, resolution=(nr, ntheta//2),
                center=(0, 0, 0), normal=(0, 0, 1),
                reflection_symmetry=False, axial_symmetry=False, name=f"half_of_{name}")
        mesh = ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=name)

    elif axial_symmetry:
        mesh_slice = mesh_disk(radius=radius, _theta_max=_theta_max/ntheta, resolution=(nr, 1),
                center=(0, 0, 0), normal=(0, 0, 1),
                reflection_symmetry=False, axial_symmetry=False, name=f"slice_of_{name}")
        mesh = AxialSymmetricMesh(mesh_slice, axis=Oz_axis, nb_repetitions=ntheta - 1, name=name)

    else:
        theta_range = np.linspace(0, _theta_max, ntheta+1)
        r_range = np.linspace(0.0, radius, nr+1)
        nodes = np.array([(0, r*sin(t), r*cos(t)) for (r, t) in product(r_range, theta_range)])
        panels = np.array([(j+i*(ntheta+1), j+1+i*(ntheta+1), j+1+(i+1)*(ntheta+1), j+(i+1)*(ntheta+1))
                                for (i, j) in product(range(0, nr), range(0, ntheta))])

        mesh = Mesh(nodes, panels, name=name)

    mesh.heal_mesh()
    mesh.translate(center)
    mesh.rotate_around_center_to_align_vectors(center, mesh.faces_normals[0], normal)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh


def mesh_vertical_cylinder(*, length=10.0, radius=1.0, center=(0, 0, 0),
        resolution=(2, 8, 10), axial_symmetry=False, reflection_symmetry=False, name=None, _theta_max=2*pi):
    """Vertical cylinder.

    Total number of panels = (2*resolution[0] + resolution[2])*resolution[1]

    Parameters
    ----------
    length : float, optional
        length of the cylinder
    radius : float, optional
        radius of the cylinder
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the cylinder
    resolution : 3-ple of int, optional
        (number of panel along a radius at the end, number of panels around a slice, number of slices)
        Mnemonic: same ordering as the cylindrical coordinates (nr, ntheta, nz)
    axial_symmetry : bool, optional
        if True, returns an AxialSymmetricMesh
    reflection_symmetry : bool, optional
        if True, returns a ReflectionSymmetricMesh
    name : str, optional
        a string naming the mesh
    _theta_max: float, optional
        internal parameter, to return an arc circle instead of a full circle
    """
    assert length > 0, "Length of a cylinder should be given as a positive value."
    assert radius > 0, "Radius of a cylinder should be given as a positive value."

    assert len(resolution) == 3, "Resolution of a cylinder should be given as a 3-ple of values."
    assert all([h >= 0 for h in resolution]), "Resolution of a cylinder should be given as positive values."
    assert all([i == int(i) for i in resolution]), "Resolution of a cylinder should be given as integer values."

    assert len(center) == 3, "Position of the center of a cylinder should be given a 3-ple of values."

    nr, ntheta, nz = resolution

    if name is None:
        name = f"cylinder_{next(Mesh._ids)}"

    LOG.debug(f"New vertical cylinder of length {length}, radius {radius} and resolution {resolution}, named {name}.")

    if reflection_symmetry and axial_symmetry:
        raise NotImplementedError("Vertical cylinders with both symmetries have not been implemented.")

    if reflection_symmetry:
        if ntheta % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        half_cylinder = mesh_vertical_cylinder(length=length, radius=radius, center=(0, 0, 0),
                resolution=(nr, ntheta//2, nz), reflection_symmetry=False, axial_symmetry=False,
                name=f"half_{name}", _theta_max=_theta_max/2)

        mesh = ReflectionSymmetricMesh(half_cylinder, plane=xOz_Plane, name=name)

    elif axial_symmetry:

        mesh_slice = mesh_vertical_cylinder(length=length, radius=radius, resolution=(nr, 1, nz), center=(0, 0, 0),
                reflection_symmetry=False, axial_symmetry=False, name=f"slice_of_{name}", _theta_max=_theta_max/ntheta)
        mesh = AxialSymmetricMesh(mesh_slice, axis=Oz_axis, nb_repetitions=ntheta - 1, name=name)

    else:
        theta_range = np.linspace(0, _theta_max, ntheta+1)
        z_range = np.linspace(-length/2, length/2, nz+1)
        if nr > 0:
            r_range = np.linspace(0.0, radius, nr+1)
            nodes = np.concatenate([
                np.array([(r*sin(t), r*cos(t), -length/2) for (r, t) in product(r_range, theta_range)]),
                np.array([(radius*sin(t), radius*cos(t), z) for (z, t) in product(z_range, theta_range)]),
                np.array([(r*sin(t), r*cos(t), length/2) for (r, t) in product(r_range[::-1], theta_range)]),
                ])
        else:
            r_range = np.array([])
            nodes = np.array([(radius*sin(t), radius*cos(t), z) for (z, t) in product(z_range, theta_range)])
        panels = np.array([(j+i*(ntheta+1), j+(i+1)*(ntheta+1), j+1+(i+1)*(ntheta+1), j+1+i*(ntheta+1), )
                                for (i, j) in product(range(nz+2*(len(r_range))), range(ntheta))])

        mesh = Mesh(nodes, panels, name=name)

    mesh.heal_mesh()
    mesh.translate(center)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh


def mesh_horizontal_cylinder(*, length=10.0, radius=1.0, center=(0, 0, 0),
        resolution=(2, 8, 10), reflection_symmetry=False, translation_symmetry=False, name=None, _theta_max=2*pi):
    """Cylinder aligned along Ox axis.

    Total number of panels = (2*resolution[0] + resolution[2])*resolution[1]

    Parameters
    ----------
    length : float, optional
        length of the cylinder
    radius : float, optional
        radius of the cylinder
    center : 3-ple or array of shape (3,), optional
        position of the geometric center of the cylinder
    resolution : 3-ple of int, optional
        (number of panel along a radius at the end, number of panels around a slice, number of slices)
        Mnemonic: same ordering as the cylindrical coordinates (nr, ntheta, nz)
    reflection_symmetry : bool, optional
        if True, returns a ReflectionSymmetricMesh
    translation_symmetry : bool, optional
        if True, uses a TranslationalSymmetricMesh internally for the main part of the cylinder
    name : str, optional
        a string naming the mesh
    _theta_max: float, optional
        internal parameter, to return an arc circle instead of a full circle
    """

    assert length > 0, "Length of a cylinder should be given as a positive value."
    assert radius > 0, "Radius of a cylinder should be given as a positive value."

    assert len(resolution) == 3, "Resolution of a cylinder should be given as a 3-ple of values."
    assert all([h >= 0 for h in resolution]), "Resolution of a cylinder should be given as positive values."
    assert all([i == int(i) for i in resolution]), "Resolution of a cylinder should be given as integer values."

    assert len(center) == 3, "Position of the center of a cylinder should be given a 3-ple of values."

    if name is None:
        name = f"cylinder_{next(Mesh._ids)}"

    LOG.debug(f"New horizontal cylinder of length {length}, radius {radius} and resolution {resolution}, named {name}.")

    nr, ntheta, nx = resolution

    if reflection_symmetry:
        if ntheta % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        half_cylinder = mesh_horizontal_cylinder(
                length=length, radius=radius, center=(0, 0, 0),
                resolution=(nr, ntheta//2, nx),
                reflection_symmetry=False, translation_symmetry=translation_symmetry,
                name=f"half_{name}", _theta_max=_theta_max/2,
                )

        mesh = ReflectionSymmetricMesh(half_cylinder, plane=xOz_Plane, name=name)

    else:
        if translation_symmetry:
            slice = mesh_horizontal_cylinder(
                    length=length/nx, radius=radius, center=(-length/2 + length/(2*nx), 0, 0),
                    resolution=(0, ntheta, 1),
                    reflection_symmetry=False, translation_symmetry=False,
                    name=f"slice_of_{name}", _theta_max=_theta_max,
                    )

            open_cylinder = TranslationalSymmetricMesh(
                    slice, translation=np.asarray([length / nx, 0.0, 0.0]),
                    nb_repetitions=nx-1, name=f"open_{name}")

        else: # General case
            theta_range = np.linspace(0, _theta_max, ntheta+1)
            x_range = np.linspace(-length/2, length/2, nx+1)
            nodes = np.array([(x, radius*sin(t), -radius*cos(t)) for (x, t) in product(x_range, theta_range)])

            panels = np.array([(i+j*(ntheta+1), i+1+j*(ntheta+1), i+1+(j+1)*(ntheta+1), i+(j+1)*(ntheta+1))
                                for (i, j) in product(range(ntheta), range(nx))])

            open_cylinder = Mesh(nodes, panels, name=f"open_{name}")

        if nr > 0:
            side = mesh_disk(radius=radius, center=(-length/2, 0, 0), normal=(-1, 0, 0),
                             reflection_symmetry=False, resolution=(nr, ntheta), name=f"side_of_{name}",
                             _theta_max=_theta_max)
            other_side = side.mirrored(yOz_Plane, name=f"other_side_of_{name}")
            mesh = CollectionOfMeshes([open_cylinder, side, other_side], name=name)
            if not translation_symmetry:
                mesh = mesh.merged()
        else:
            mesh = open_cylinder.copy(name=name)

    mesh.heal_mesh()
    mesh.translate(center)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh

