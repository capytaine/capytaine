#!/usr/bin/env python
# coding: utf-8
"""Tests for the mesh import of mesh from external libraries"""

import os
import pytest
import logging
from unittest import mock
from pathlib import Path

import numpy as np
import xarray as xr

from capytaine.io.mesh_writers import write_STL, write_GDF
from capytaine.io.mesh_loaders import load_STL, load_HST, load_GDF
import capytaine as cpt

try:
    import meshio
    import pygmsh
    import gmsh
except ImportError:
    pygmsh = None

try:
    import vtk
except ImportError:
    vtk = None

try:
    import h5py
except ImportError:
    h5py = None

offset=1e-2 # small offset so you can clip at z=0
omega = 0.5*2*np.pi
rho_water = 1e3

@pytest.mark.skipif(vtk is None,
                    reason='vtk is not installed')
def test_STL(tmp_path):
    mesh = cpt.Sphere().mesh.merged()
    filepath = tmp_path / "test.stl"
    write_STL(filepath, mesh.vertices, mesh.faces)
    reloaded_mesh = load_STL(str(filepath), name="Bla")

    assert reloaded_mesh.name == "Bla"
    assert np.allclose(mesh.vertices, reloaded_mesh.vertices)

#################################################################

def parse_hst(file_content):
    # Workaround to avoid actually writing/reading a file on disk
    with mock.patch('capytaine.io.mesh_loaders._check_file', lambda foo: None):
        with mock.patch('capytaine.io.mesh_loaders.open', mock.mock_open(read_data=file_content)):
            return load_HST("mocked/filename")


def test_hst_type_0_triangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_type_0_quadrangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2/3, 2/3, 0.0]]))


def test_hst_type_1_quadrangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_type_1_triangle():
    mesh = parse_hst("""
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 1
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[2/3, 2/3, 0.0]]))


def test_hst_implicit_coordinate_numbering():
    mesh = parse_hst("""
        COORDINATES
        0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
        """)
    assert mesh.vertices.shape == (4, 3)
    assert mesh.faces.shape == (1, 4)
    assert np.allclose(mesh.faces_centers, np.array([[1/2, 1/2, 0.0]]))


def test_hst_coordinate_numbering_error():
    with pytest.raises(ValueError):
        mesh = parse_hst("""
            COORDINATES
            1 0.0 0.0 0.0
            2 1.0 0.0 0.0
            3 1.0 1.0 0.0
            5 0.0 1.0 0.0
            ENDCOORDINATES
            """)

def test_hst_symmetry_1():
    mesh = parse_hst("""
        SYMMETRY 1
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.Mesh)
    assert mesh.merged().nb_faces == 2


def test_hst_symmetry_2():
    mesh = parse_hst("""
        SYMMETRY 2
        COORDINATES
        1 0.0 0.0 0.0
        2 1.0 0.0 0.0
        3 1.0 1.0 0.0
        4 0.0 1.0 0.0
        ENDCOORDINATES
        PANEL TYPE 0
        1 2 3 4
        ENDPANEL
            """)
    assert isinstance(mesh, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half, cpt.ReflectionSymmetricMesh)
    assert isinstance(mesh.half.half, cpt.Mesh)
    assert mesh.merged().nb_faces == 4


def test_hst_ignored_lines(caplog):
    with caplog.at_level(logging.WARNING):
        mesh = parse_hst("""
            USER mancellin
            NBBODY 1
            COORDINATES
            1 0.0 0.0 0.0
            2 1.0 0.0 0.0
            3 1.0 1.0 0.0
            4 0.0 1.0 0.0
            ENDCOORDINATES
            PANEL TYPE 0
            1 2 3 4
            ENDPANEL
                """)
    assert "HST mesh reader ignored" in caplog.text
    assert "USER mancellin" in caplog.text
    assert "NBBODY 1" in caplog.text

#################################################################

def generate_wavebot():
    T1=0.16
    T2=0.37
    r1=0.88
    r2=0.35
    with pygmsh.occ.Geometry() as geom:
        cyl = geom.add_cylinder([0,0,0], [0,0,-T1], r1)
        cone = geom.add_cone([0,0,-T1], [0,0,-T2], r1, r2)
        geom.translate(cyl, [0, 0, offset])
        geom.translate(cone, [0, 0, offset])
        geom.boolean_union([cyl, cone])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = np.pi*T1*r1**2 + 1/3 * np.pi * T2 * (r2**2 + r2 * r1 + r1**2)
    return (mesh, vol_exp)

def generate_cylinder():
    T=0.52
    r1=0.88
    with pygmsh.occ.Geometry() as geom:
        cyl = geom.add_cylinder([0,0,0], [0,0,-T], r1)
        geom.translate(cyl, [0, 0, offset])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = np.pi*r1**2*T
    return (mesh, vol_exp)

def generate_sphere():
    r1=0.88
    with pygmsh.occ.Geometry() as geom:
        sphere = geom.add_ball(center=[0,0,0], radius=r1)
        geom.translate(sphere, [0, 0, offset])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = 1/2*np.pi*4/3*r1**3
    return (mesh, vol_exp)


@pytest.mark.skipif(vtk is None or pygmsh is None,
                    reason='Neither vtk nor meshio are installed')
def test_write_and_read_STL(tmp_path):
    mesh, _ = generate_cylinder()
    from capytaine.io.meshio import load_from_meshio
    cpt_mesh = load_from_meshio(mesh)
    cpt_mesh.heal_mesh()

    # Write with Meshio and reload with Meshmagick
    mesh.write(tmp_path / "wavebot.stl")
    cpt_mesh_2 = cpt.load_mesh(tmp_path / "wavebot.stl")
    assert np.allclose(cpt_mesh.vertices, cpt_mesh_2.vertices)

    # # Write with Meshmagick and reload with Meshio
    # # FAILING FOR NOW
    # # ERROR IN STL LOADER IN MESHMAGICK
    # write_STL(tmp_path / "wavebot2.stl", fb.mesh.vertices, fb.mesh.faces)
    # mesh2 = meshio.read(tmp_path / "wavebot2.stl")
    # fb3 = cpt.FloatingBody.from_meshio(mesh2)
    # assert np.allclose(fb.mesh.vertices, fb3.mesh.vertices)


@pytest.mark.skipif(pygmsh is None,
                    reason='pygmsh and/or meshio is not installed')
@pytest.mark.parametrize("generate_pygmsh",
                         [generate_cylinder, generate_sphere, generate_wavebot])
def test_from_meshio_pygmsh(generate_pygmsh, tmp_path):
    mesh, vol_exp = generate_pygmsh()
    fb = cpt.FloatingBody(mesh=mesh)
    fb.mesh.heal_mesh()
    fb.keep_immersed_part()

    vol = fb.mesh.volume
    assert pytest.approx(vol_exp, rel=1e-1) == vol

    fb.add_translation_dof((0,0,1))

    test_matrix = xr.Dataset(coords={
        'rho': rho_water,
        'water_depth': [np.infty],
        'omega': omega,
        'wave_direction': 0,
        'radiating_dof': list(fb.dofs.keys()),
        })

    solver = cpt.BEMSolver()
    data = solver.fill_dataset(test_matrix,
                               bodies=[fb],
                               hydrostatics=True,
                               mesh=True,
                               wavelength=True,
                               wavenumber=True)


def test_write_and_load_gdf(tmpdir):
    mesh_path = tmpdir.join("temp_mesh.gdf")

    original_mesh = cpt.mesh_horizontal_cylinder()
    write_GDF(str(mesh_path), original_mesh.vertices, original_mesh.faces, ulen=1, gravity=9.81, isx=0, isy=0)

    read_mesh = load_GDF(str(mesh_path))

    np.testing.assert_allclose(
        read_mesh.vertices[read_mesh.faces],
        original_mesh.vertices[original_mesh.faces],
        atol=1e-6
        )

    np.testing.assert_allclose(
        read_mesh.faces_areas,
        original_mesh.faces_areas
        )

    np.testing.assert_allclose(
        read_mesh.faces_normals,
        original_mesh.faces_normals,
        atol=1e-6
        )

    np.testing.assert_allclose(
        read_mesh.faces_centers,
        original_mesh.faces_centers
        )

    np.testing.assert_allclose(
        read_mesh.faces_radiuses,
        original_mesh.faces_radiuses
        )

    np.testing.assert_allclose(
        read_mesh.volume,
        original_mesh.volume
        )


@pytest.mark.skipif(h5py is None,
                    reason='h5py and/or meshio is not installed')
def test_MED_file():
    mesh = cpt.load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/barge.med"))
    assert mesh.nb_faces == 187


def test_pnl(tmpdir):
    mesh_path = tmpdir.join("temp_mesh.pnl")
    mesh = cpt.mesh_sphere()
    from capytaine.io.mesh_writers import write_PNL
    write_PNL(mesh_path, mesh.vertices, mesh.faces)
    reloaded_mesh = cpt.load_mesh(mesh_path)
    np.testing.assert_equal(mesh.faces, reloaded_mesh.faces)
    np.testing.assert_allclose(mesh.vertices, reloaded_mesh.vertices, atol=1e-5)
