# Copyright 2025 Mews Labs
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

from functools import lru_cache
from pathlib import Path

import numpy as np
import xarray as xr

from capytaine.new_meshes import Mesh
from capytaine.new_meshes.io import load_mesh

@lru_cache
def simple_triangle_mesh():
    """Simple triangle mesh for testing."""
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    faces = np.array([[0, 1, 2]], dtype=int)
    return Mesh(vertices, faces)


@lru_cache
def simple_quad_mesh():
    """Simple quad mesh for testing."""
    vertices = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float
    )
    faces = np.array([[0, 1, 2, 3]], dtype=int)
    return Mesh(vertices, faces)


def test_mesh_to_xarray_triangle():
    mesh = simple_triangle_mesh()
    ds = mesh.export_to_xarray()
    assert isinstance(ds, xr.Dataset)
    assert "mesh_vertices" in ds.data_vars
    assert ds["mesh_vertices"].shape == (1, 3, 3)


def test_mesh_to_xarray_quad():
    mesh = simple_quad_mesh()
    ds = mesh.export_to_xarray()
    assert isinstance(ds, xr.Dataset)
    assert "mesh_vertices" in ds.data_vars
    assert ds["mesh_vertices"].shape == (1, 4, 3)


def test_mesh_to_xarray_roundtrip():
    mesh = simple_triangle_mesh()
    ds = mesh.export_to_xarray()
    mesh_ = load_mesh(ds)
    assert np.allclose(mesh_.faces_centers, mesh.faces_centers)
    assert np.allclose(mesh_.faces_normals, mesh.faces_normals)


def test_xarray_to_mesh_from_file_detecting_extension(tmpdir):
    mesh = simple_triangle_mesh()
    ds = mesh.export_to_xarray()
    filepath = Path(tmpdir / "mesh.nc")
    ds.to_netcdf(filepath)
    mesh_ = load_mesh(filepath)
    assert np.allclose(mesh_.faces_centers, mesh.faces_centers)
    assert np.allclose(mesh_.faces_normals, mesh.faces_normals)


def test_xarray_to_mesh_from_file_explicit_format(tmpdir):
    mesh = simple_triangle_mesh()
    ds = mesh.export_to_xarray()
    filepath = Path(tmpdir / "mesh")
    ds.to_netcdf(filepath)
    mesh_ = load_mesh(filepath, file_format="nc")
    assert np.allclose(mesh_.faces_centers, mesh.faces_centers)
    assert np.allclose(mesh_.faces_normals, mesh.faces_normals)
