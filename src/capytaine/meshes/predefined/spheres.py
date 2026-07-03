# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Generate spherical bodies."""

import logging

import numpy as np
from numpy import pi

from capytaine.meshes.symmetric_meshes import RotationSymmetricMesh

LOG = logging.getLogger(__name__)


def mesh_sphere(*, radius=1.0, center=(0.0, 0.0, 0.0),
        resolution=(10, 10), faces_max_radius=None,
        axial_symmetry=False, name=None):
    """Sphere

    Parameters
    ----------
    radius : float
        radius of the sphere
    center : 3-ple or array of shape (3,)
        position of the geometric center of the sphere
    resolution : couple of ints
        number of panels along a meridian (or number of parallels-1) and
        along a parallel (or number of meridians-1)
    faces_max_radius : float, optional
        maximal radius of a panel. (Default: no maximal radius.)
        If the provided resolution is too coarse, the number of panels is
        changed to fit the constraint on the maximal radius.
    axial_symmetry : bool
        if True, use the axial symmetry to build the mesh (default: False)
    name : string
        a name identifying the sphere (default: "sphere_id" where id is an unique integer).
    """

    ntheta, nphi = resolution
    if faces_max_radius is not None:
        perimeter = 2*np.pi*radius
        estimated_max_radius = np.hypot(perimeter/ntheta, perimeter/nphi)/2
        if estimated_max_radius > faces_max_radius:
            ntheta = nphi = int(np.ceil(perimeter / (np.sqrt(2)*faces_max_radius)))

    theta = np.linspace(0.0, pi, ntheta+1)
    points_on_a_meridian = radius * np.stack([np.sin(theta), np.zeros_like(theta), -np.cos(theta)], axis=1)

    mesh = RotationSymmetricMesh.from_profile_points(points_on_a_meridian, n=nphi, name=name)

    if not axial_symmetry:
        mesh = mesh.merged()

    mesh = mesh.translated(center, name=name)
    return mesh
