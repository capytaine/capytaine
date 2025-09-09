"""Generate spherical bodies."""
# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

import numpy as np

from capytaine.meshes import Mesh
from capytaine.meshes.predefined import mesh_sphere
from capytaine.bodies.bodies import FloatingBody

LOG = logging.getLogger(__name__)


class Sphere(FloatingBody):
    """Sphere
    Deprecated: please prefer capytaine.meshes.predefined.mesh_sphere()

    Parameters
    ----------
    radius : float
        radius of the sphere
    center : 3-ple or array of shape (3,)
        position of the geometric center of the sphere
    ntheta : int
        number of panels along a meridian (or number of parallels-1)
    nphi : int
        number of panels along a parallel (or number of meridians-1)
    axial_symmetry : bool
        if True, use the axial symmetry to build the mesh (default: True)
    clip_free_surface : bool
        if True, only mesh the part of the sphere where z < 0 (default: False),
        can be used with center to obtain any clipped sphere,
        if True, then ntheta is the number of parallel below the free surface.
    name : string
        a name identifying the sphere (default: "sphere_id" where id is an unique integer).
    """

    def __init__(self, *, radius=1.0, center=(0, 0, 0),
                 ntheta=10, nphi=10, clip_free_surface=False,
                 axial_symmetry=True, clever=None,
                 name=None):

        LOG.warning("Deprecation warning: The class Sphere() is deprecated. "
                "Please prefer the function capytaine.meshes.predefined.mesh_sphere()")

        if clever is not None:
            LOG.warning("Deprecation warning: `clever` argument for Sphere is deprecated. "
                        "Use `axial_symmetry` instead.")

        if name is None:
            name = f"sphere_{next(Mesh._ids)}"

        if clip_free_surface:
            if center[2] < -radius:  # fully immersed
                pass
            elif center[2] < radius:
                ntheta = int(ntheta*np.pi/np.arccos(center[2]/radius))
            else:
                raise ValueError("Impossible to mesh the immersed hull of a sphere completely out of the water")

        mesh = mesh_sphere(radius=radius, center=center, resolution=(ntheta, nphi), axial_symmetry=axial_symmetry, name=f"{name}_mesh")

        if clip_free_surface:
            mesh.keep_immersed_part()

        self.radius = radius
        self.geometric_center = np.array(center, dtype=float)
        FloatingBody.__init__(self, mesh=mesh, name=name)
