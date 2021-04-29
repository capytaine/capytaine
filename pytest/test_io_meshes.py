#!/usr/bin/env python
# coding: utf-8
"""Tests for the mesh submodule: definition and transformation of base meshes."""

import pytest

import numpy as np

from capytaine.io.mesh_writers import write_STL
from capytaine.io.mesh_loaders import load_STL
from capytaine.bodies.predefined import Sphere

try:
    import vtk
except ImportError:
    vtk = None

@pytest.mark.skipif(vtk is None,
                    reason='vtk is not installed')
def test_STL(tmp_path):
    mesh = Sphere().mesh.merged()
    filepath = tmp_path / "test.stl"
    write_STL(filepath, mesh.vertices, mesh.faces)
    reloaded_mesh = load_STL(str(filepath), name="Bla")

    assert reloaded_mesh.name == "Bla"
    assert np.allclose(mesh.vertices, reloaded_mesh.vertices)
    # Cannot compare the faces. The STL writer changed all quadrangles to two triangles.
