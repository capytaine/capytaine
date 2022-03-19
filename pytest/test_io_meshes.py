#!/usr/bin/env python
# coding: utf-8
"""Tests for the mesh import of mesh from external libraries"""

import pytest

import numpy as np
import xarray as xr

from capytaine.io.mesh_writers import write_STL
from capytaine.io.mesh_loaders import load_STL
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
    fb = cpt.FloatingBody.from_meshio(mesh)
    fb.mesh.heal_mesh()

    # Write with Meshio and reload with Meshmagick
    mesh.write(tmp_path / "wavebot.stl")
    fb2 = cpt.FloatingBody.from_file(str(tmp_path / "wavebot.stl"))
    assert np.allclose(fb.mesh.vertices, fb2.mesh.vertices)

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
    fb = cpt.FloatingBody.from_meshio(mesh)
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

