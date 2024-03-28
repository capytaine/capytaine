"""Generate rectangular bodies."""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
import logging
from itertools import product

import numpy as np

from capytaine.meshes.geometry import xOz_Plane, yOz_Plane
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric import TranslationalSymmetricMesh, ReflectionSymmetricMesh
from capytaine.meshes.collections import CollectionOfMeshes

LOG = logging.getLogger(__name__)

def mesh_rectangle(*, size=(5.0, 5.0), center=(0.0, 0.0, 0.0),
        resolution=(5, 5), faces_max_radius=None,
        normal=(0.0, 0.0, 1.0),
        translation_symmetry=False, reflection_symmetry=False,
        name=None):
    """One-sided rectangle.

    By default, the rectangle is horizontal, the normals are oriented upwards.

    Parameters
    ----------
    size : couple of floats, optional
        dimensions of the rectangle (width and height)
    center : 3-ple of floats, optional
        position of the geometric center of the rectangle, default: (0, 0, 0)
    resolution : couple of ints, optional
        number of faces along each of the two directions
    faces_max_radius : float, optional
        maximal radius of a panel. (Default: no maximal radius.)
        If the provided resolution is too coarse, the number of panels is
        changed to fit the constraint on the maximal radius.
    normal: 3-ple of floats, optional
        normal vector, default: (0, 0, 1)
    translation_symmetry : bool, optional
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

    width, height = size
    nw, nh = resolution

    if faces_max_radius is not None:
        estimated_max_radius = np.hypot(width/nw, height/nh)/2
        if estimated_max_radius > faces_max_radius:
            nw = int(np.ceil(width / (np.sqrt(2) * faces_max_radius)))
            nh = int(np.ceil(height / (np.sqrt(2) * faces_max_radius)))

    if name is None:
        name = f"rectangle_{next(Mesh._ids)}"

    if translation_symmetry and reflection_symmetry:
        raise NotImplementedError("Rectangle generation with both reflection and translation symmetries "
                                  "has not been implemented.")

    if reflection_symmetry:
        if nw % 2 == 1:
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        half_mesh = mesh_rectangle(size=(width/2, height), resolution=(nw//2, nh),
                                   center=(0, -width/4, 0), normal=(0.0, 0.0, 1.0),
                                   translation_symmetry=False, reflection_symmetry=False,
                                   name=f"half_of_{name}")
        mesh = ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=name)

    elif translation_symmetry:
            strip = mesh_rectangle(size=(width/nw, height), resolution=(1, nh),
                                   center=(0, -width/2 + width/(2*nw), 0), normal=(0.0, 0.0, 1.0),
                                   translation_symmetry=False, reflection_symmetry=False,
                                   name=f"strip_of_{name}")
            mesh = TranslationalSymmetricMesh(strip,
                                              translation=np.asarray([0, width/nw, 0]),
                                              nb_repetitions=int(nw)-1,
                                              name=name)

    else:
        y_range = np.linspace(-width/2, width/2, nw+1)
        z_range = np.linspace(-height/2, height/2, nh+1)
        nodes = np.array(list(product([0.0], y_range, z_range)), dtype=float)
        panels = np.array([(j+i*(nh+1), j+1+i*(nh+1), j+1+(i+1)*(nh+1), j+(i+1)*(nh+1))
                             for (i, j) in product(range(nw), range(nh))])

        mesh = Mesh(nodes, panels, name=name)

    mesh.heal_mesh()
    mesh.translate(center)
    mesh.rotate_around_center_to_align_vectors(center, mesh.faces_normals[0], normal)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh


def mesh_parallelepiped(size=(1.0, 1.0, 1.0), center=(0, 0, 0),
                        resolution=(4, 4, 4), faces_max_radius=None,
                        missing_sides=set(), reflection_symmetry=False, translation_symmetry=False,
                        name=None):
    """Six rectangles forming a parallelepiped.

    Parameters
    ----------
    size : 3-ple of floats, optional
        dimensions of the parallelepiped (width, thickness, height) for coordinates (x, y, z).
    center : 3-ple of floats, optional
        coordinates of the geometric center of the parallelepiped
    resolution : 3-ple of ints, optional
        number of faces along the three directions
    faces_max_radius : float, optional
        maximal radius of a panel. (Default: no maximal radius.)
        If the provided resolution is too coarse, the number of panels is
        changed to fit the constraint on the maximal radius.
    missing_sides : set of string, optional
        if one of the keyword "top", "bottom", "front", "back", "left", "right" is in the set,
        then the corresponding side is not included in the parallelepiped.
        May be ignored when building a mesh with a symmetry.
    reflection_symmetry : bool, optional
        use xOz and yOz symmetry plane to generate the mesh
    translation_symmetry : bool, optional
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

    width, thickness, height = size
    nw, nth, nh = resolution

    if faces_max_radius is not None:
        dw, dh, dth = width/nw, height/nh, thickness/nth
        estimated_max_radius = max(
                np.hypot(dw, dh)/2,
                np.hypot(dw, dth)/2,
                np.hypot(dth, dh)/2,
                )
        if estimated_max_radius > faces_max_radius:
            nw = int(np.ceil(width / (np.sqrt(2) * faces_max_radius)))
            nth = int(np.ceil(thickness / (np.sqrt(2) * faces_max_radius)))
            nh = int(np.ceil(height / (np.sqrt(2) * faces_max_radius)))

    if name is None:
        name = f"rectangular_parallelepiped_{next(Mesh._ids)}"

    if translation_symmetry and reflection_symmetry:
        raise NotImplementedError("Parallelepiped generation with both reflection and translation symmetries "
                                  "has not been implemented.")

    if reflection_symmetry:
        if (nw % 2 == 1 or nth % 2 == 1):
            raise ValueError("To use the reflection symmetry of the mesh, "
                             "it should have an even number of panels in this direction.")

        missing_sides_in_quarter = missing_sides | {"right", "back"}
        quarter_mesh = mesh_parallelepiped(
                size=(width/2, thickness/2, height), resolution=(nw//2, nth//2, nh),
                center=(-width/4, -thickness/4, 0), missing_sides=missing_sides_in_quarter,
                reflection_symmetry=False, translation_symmetry=False,
                name=f"quarter_of_{name}"
                )

        half_mesh = ReflectionSymmetricMesh(quarter_mesh, plane=yOz_Plane, name=f"half_of_{name}")
        mesh = ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=f"{name}")

    elif translation_symmetry:

        missing_sides_in_strip = missing_sides | {"left", "right"}
        strip = mesh_parallelepiped(
                size=(width/nw, thickness, height), resolution=(1, nth, nh),
                center=(-width/2 + width/(2*nw), 0, 0), missing_sides=missing_sides_in_strip,
                reflection_symmetry=False, translation_symmetry=False,
                name=f"strip_of_{name}"
                )

        open_parallelepiped = TranslationalSymmetricMesh(
            strip,
            translation=(width/nw, 0, 0), nb_repetitions=int(nw)-1,
            name=f"body_of_{name}"
        )

        components_of_mesh = [open_parallelepiped]
        if "right" not in missing_sides:
            components_of_mesh.append(
                    mesh_rectangle(
                        size=(thickness, height), resolution=(nth, nh),
                        center=(width/2, 0, 0), normal=(1, 0, 0),
                        name=f"right_side_of_{name}"
                        ))
        if "left" not in missing_sides:
            components_of_mesh.append(
                    mesh_rectangle(
                        size=(thickness, height), resolution=(nth, nh),
                        center=(-width/2, 0, 0), normal=(-1, 0, 0),
                        name=f"left_side_of_{name}"
                        ))

        mesh = CollectionOfMeshes(components_of_mesh, name=name)

    else:

        sides = []
        if "left" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(thickness, height), resolution=(nth, nh), center=(-width/2, 0, 0),
                        normal=(-1, 0, 0), name=f"left_of_{name}")
                    )
        if "right" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(thickness, height), resolution=(nth, nh), center=(width/2, 0, 0),
                        normal=(1, 0, 0), name=f"right_of_{name}")
                    )
        if "front" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(width, height), resolution=(nw, nh), center=(0, -thickness/2, 0),
                        normal=(0, -1, 0), name=f"front_of_{name}")
                    )
        if "back" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(width, height), resolution=(nw, nh), center=(0, thickness/2, 0),
                        normal=(0, 1, 0), name=f"back_of_{name}")
                    )
        if "top" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(thickness, width), resolution=(nth, nw), center=(0, 0, height/2),
                        normal=(0, 0, 1), name=f"top_of_{name}")
                    )
        if "bottom" not in missing_sides:
            sides.append(
                    mesh_rectangle(size=(thickness, width), resolution=(nth, nw), center=(0, 0, -height/2),
                        normal=(0, 0, -1), name=f"bottom_of_{name}")
                    )
        mesh = CollectionOfMeshes(sides, name=name).merged()

    mesh.heal_mesh()
    mesh.translate(center)
    mesh.geometric_center = np.asarray(center, dtype=float)
    return mesh
