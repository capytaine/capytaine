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
"""Helper functions shared by several test files.

Meshes and Green functions are cached, since they are not modified by the tests.
Bodies and solvers are built anew at each call, since they are mutable (dofs, timer, cached matrices, ...).
"""

from functools import lru_cache

import numpy as np

import capytaine as cpt
from capytaine.meshes import ReflectionSymmetricMesh, RotationSymmetricMesh


################################
#  Green functions and solvers #
################################

@lru_cache
def green_function(**kwargs):
    """Shared Delhommeau instance, to avoid reloading the tabulation in each test."""
    return cpt.Delhommeau(**kwargs)


def solver(method="indirect", **gf_kwargs):
    """New BEMSolver (with its own timer and cached matrices) using a shared Green function."""
    return cpt.BEMSolver(green_function=green_function(**gf_kwargs), method=method)


class BrokenGreenFunction:
    """Green function failing for wavenumber < 2.0, to test the handling of failed resolutions."""
    floating_point_precision = None
    exportable_settings = {}

    def evaluate(self, m1, m2, *, wavenumber, **kwargs):
        if wavenumber < 2.0:
            raise NotImplementedError("I'm potato")
        else:
            return green_function().evaluate(m1, m2, wavenumber=wavenumber, **kwargs)


def broken_bem_solver():
    return cpt.BEMSolver(green_function=BrokenGreenFunction())


############
#  Meshes  #
############

@lru_cache
def small_sphere_mesh():
    """Coarse immersed half sphere of radius 1 centered on the free surface."""
    return cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()


def small_sphere_body(dofs=("Surge",), *, lid=False, name=None):
    """New FloatingBody with the mesh `small_sphere_mesh()`.

    Parameters
    ----------
    dofs: "rigid" or iterable of str
        Either all six rigid body dofs (around the origin) or the names of the subset of them.
    lid: bool
        Add a lid on the free surface for irregular frequency removal.
    name: str, optional
    """
    mesh = small_sphere_mesh()
    lid_mesh = mesh.generate_lid() if lid else None
    if dofs == "rigid":
        dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    else:
        dofs = cpt.rigid_body_dofs(only=dofs, rotation_center=(0, 0, 0))
    return cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=dofs, center_of_mass=(0, 0, 0), name=name)


@lru_cache
def single_panel():
    """A single non-planar panel, to build small meshes with symmetries."""
    vertices = np.array([[0.5, 0.0, 0.0], [0.5, 0.0, -0.5], [0.5, 0.5, -0.3], [0.5, 0.5, -0.2]])
    return cpt.Mesh(vertices=vertices, faces=np.array([[0, 1, 2, 3]]))


@lru_cache
def symmetric_meshes_of_single_panel():
    return {
        "reflection": ReflectionSymmetricMesh(single_panel(), plane="xOz"),
        "nested_reflections": ReflectionSymmetricMesh(ReflectionSymmetricMesh(single_panel(), plane="xOz"), plane="yOz"),
        "rotation_2": RotationSymmetricMesh(single_panel(), n=2, axis='z+'),
        "rotation_3": RotationSymmetricMesh(single_panel(), n=3, axis='z+'),
        "rotation_4": RotationSymmetricMesh(single_panel(), n=4, axis='z+'),
        "rotation_5": RotationSymmetricMesh(single_panel(), n=5, axis='z+'),
        "dihedral": RotationSymmetricMesh(ReflectionSymmetricMesh(single_panel(), plane="xOz"), n=3, axis='z+'),
    }
