#!/usr/bin/env python
# coding: utf-8
"""Generate meshes of cylinders and disks"""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
from itertools import product

import numpy as np

from capytaine.meshes.geometry import xOy_Plane, xOz_Plane, yOz_Plane, e_x, e_z, Oz_axis
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import TranslationalSymmetricMesh, AxialSymmetricMesh, ReflectionSymmetricMesh

LOG = logging.getLogger(__name__)


##########
#  Disk  #
##########

def mesh_disk(*, radius=1.0, resolution=(3, 5), center=(0, 0, 0), normal=(1, 0, 0),
        theta_max=2*np.pi, reflection_symmetry=False, axial_symmetry=False, name=None):
    """(One-sided) disk.

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
    theta_max: float, optional
        if you want an arc circle instead of a full circle
        default: 2π for a full circle
    axial_symmetry : bool, optional
        if True, use the axial symmetry
        default: False
    reflection_symmetry : bool, optional
        if True, use the reflection symmetry
        default: False
    name : str, optional
        a string naming the mesh
    """
    assert radius > 0, "Radius of the disk mesh should be given as a positive value."

    assert len(resolution) == 2, "Resolution of a disk should be given as a couple a values."
    assert all([h > 0 for h in resolution]), "Resolution of the disk mesh should be given as positive values."
    assert all([i == int(i) for i in resolution]), "Resolution of a disk should be given as integer values."

    assert len(center) == 3, "Position of the center of a disk should be given a 3-ple of values."

    nr, ntheta = resolution

    if reflection_symmetry and ntheta % 2 == 1:
        raise ValueError("To use the reflection symmetry of the mesh, "
                         "it should have an even number of panels in this direction.")

    if reflection_symmetry and axial_symmetry:
        raise NotImplementedError("Disk generators with both symmetries have not been implemented.")

    if name is None:
        name = f"disk_{next(Mesh._ids)}"

    LOG.debug(f"New disk body of radius {radius} and resolution ({nr}, {ntheta}), named {name}.")

    if reflection_symmetry:
        half_mesh = mesh_disk(radius=radius, theta_max=theta_max/2, resolution=(nr, ntheta//2),
                center=(0, 0, 0), normal=(1, 0, 0),
                reflection_symmetry=False, axial_symmetry=False, name=f"half_of_{name}")
        mesh = ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=name)

    elif axial_symmetry:
        mesh_slice = mesh_disk(radius=radius, theta_max=theta_max/ntheta, resolution=(nr, 1),
                center=(0, 0, 0), normal=(1, 0, 0),
                reflection_symmetry=False, axial_symmetry=False, name=f"slice_of_{name}")
        mesh_slice.rotate_around_center_to_align_vectors((0, 0, 0), e_x, e_z)  # Convoluted way to avoid a warning message in AxialSymmetry...
        mesh = AxialSymmetricMesh(mesh_slice, axis=Oz_axis, nb_repetitions=ntheta - 1, name=name)
        mesh.rotate_around_center_to_align_vectors((0, 0, 0), e_z, e_x)

    else:
        theta_range = np.linspace(0, 2*theta_max, ntheta+1)
        r_range = np.linspace(0.0, radius, nr+1)

        nodes = np.zeros(((ntheta+1)*(nr+1), 3), dtype=float)
        for i, (r, t) in enumerate(product(r_range, theta_range)):
            y = +r * np.sin(t)
            z = -r * np.cos(t)
            nodes[i, :] = (0, y, z)

        panels = np.zeros((ntheta*nr, 4), dtype=int)

        for k, (i, j) in enumerate(product(range(0, nr), range(0, ntheta))):
            panels[k, :] = (
                j+i*(ntheta+1),
                j+1+i*(ntheta+1),
                j+1+(i+1)*(ntheta+1),
                j+(i+1)*(ntheta+1)
            )

        mesh = Mesh(nodes, panels, name=name)
        mesh.merge_duplicates()
        mesh.heal_triangles()

    mesh.rotate_around_center_to_align_vectors((0, 0, 0), mesh.faces_normals[0], normal)
    mesh.translate(center)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh


