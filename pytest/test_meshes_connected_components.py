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
import pytest
import numpy as np

import capytaine as cpt

##########################
#  Connected components  #
##########################

def test_connected_components_of_empty_mesh():
    from capytaine.meshes.geometry import connected_components
    mesh = cpt.Mesh()
    cc = connected_components(mesh)
    assert len(cc) == 0


def test_connected_components_of_sphere():
    from capytaine.meshes.geometry import connected_components
    mesh = cpt.mesh_sphere()
    cc = connected_components(mesh)
    assert len(cc) == 1


def test_connected_components_of_two_spheres():
    from capytaine.meshes.geometry import connected_components
    mesh = cpt.mesh_sphere() + cpt.mesh_sphere(center=(5, 0, 0))
    cc = connected_components(mesh)
    assert len(cc) == 2


def test_connected_components_of_open_sphere():
    from capytaine.meshes.geometry import connected_components
    mesh = cpt.mesh_sphere().immersed_part()
    cc = connected_components(mesh)
    assert len(cc) == 1


def test_connected_components_of_torus():
    from capytaine.meshes.geometry import connected_components
    pygmsh = pytest.importorskip("pygmsh")
    with pygmsh.occ.Geometry() as geom:
        geom.add_torus((0, 0, 0), 2.0, 0.5, mesh_size=0.3)
        gmsh_mesh = geom.generate_mesh(dim=2)
    mesh = cpt.load_mesh(gmsh_mesh)
    cc = connected_components(mesh)
    assert len(cc) == 1


######################################
#  Connected components at waterline #
######################################

def test_connected_components_at_waterline_of_sphere():
    from capytaine.meshes.geometry import connected_components_of_waterline
    mesh = cpt.mesh_sphere()
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 1


def test_connected_components_at_waterline_of_two_spheres():
    from capytaine.meshes.geometry import connected_components_of_waterline
    mesh = cpt.mesh_sphere() + cpt.mesh_sphere(center=(5, 0, 0))
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 2


def test_connected_components_at_waterline_of_open_sphere():
    from capytaine.meshes.geometry import connected_components_of_waterline
    mesh = cpt.mesh_sphere().immersed_part()
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 1


def test_connected_components_at_waterline_of_immersed_sphere():
    from capytaine.meshes.geometry import connected_components_of_waterline
    mesh = cpt.mesh_sphere(center=(0, 0, -5))
    cc = connected_components_of_waterline(mesh, z=0.0)
    assert len(cc) == 0


def test_connected_components_at_waterline_of_torus():
    from capytaine.meshes.geometry import connected_components_of_waterline
    pygmsh = pytest.importorskip("pygmsh")
    with pygmsh.occ.Geometry() as geom:
        geom.add_torus((0, 0, 0), 2.0, 0.5, mesh_size=0.3)
        gmsh_mesh = geom.generate_mesh(dim=2)
    mesh = cpt.load_mesh(gmsh_mesh)
    cc = connected_components_of_waterline(mesh)
    assert len(cc) == 2
