
from functools import lru_cache

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
    filepath = tmpdir / "mesh.nc"
    ds.to_netcdf(filepath)
    mesh_ = load_mesh(filepath)
    assert np.allclose(mesh_.faces_centers, mesh.faces_centers)
    assert np.allclose(mesh_.faces_normals, mesh.faces_normals)

def test_xarray_to_mesh_from_file_explicit_format(tmpdir):
    mesh = simple_triangle_mesh()
    ds = mesh.export_to_xarray()
    filepath = tmpdir / "mesh"
    ds.to_netcdf(filepath)
    mesh_ = load_mesh(filepath, file_format="nc")
    assert np.allclose(mesh_.faces_centers, mesh.faces_centers)
    assert np.allclose(mesh_.faces_normals, mesh.faces_normals)
