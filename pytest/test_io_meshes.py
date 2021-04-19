#!/usr/bin/env python
# coding: utf-8
"""Tests for the mesh submodule: definition and transformation of base meshes."""

from warnings import warn
import pytest

import numpy as np
import xarray as xr
import pygmsh

from capytaine.io.mesh_writers import write_STL
from capytaine.io.mesh_loaders import load_STL
from capytaine.bodies.predefined import Sphere
from capytaine import FloatingBody
from capytaine import BEMSolver

offset=1e-2 # small offset so you can clip at z=0

omega = 0.5*2*np.pi
rho_water = 1e3


def test_STL(tmp_path):
    try:

        mesh = Sphere().mesh.merged()
        filepath = tmp_path / "test.stl"
        write_STL(filepath, mesh.vertices, mesh.faces)
        reloaded_mesh = load_STL(str(filepath), name="Bla")

        assert reloaded_mesh.name == "Bla"
        assert np.allclose(mesh.vertices, reloaded_mesh.vertices)
        # Cannot compare the faces. The STL writer changed all quadrangles to two triangles.

    except ImportError:
        warn("VTK is not installed and thus has not been tested.")

def test_from_meshio_pygmsh_WaveBot():
    try:
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
        fb = FloatingBody.from_meshio(mesh)

        vol = fb.mesh.volume
        vol_exp = np.pi*T1*r1**2 + 1/3 * np.pi * T2 * (r2**2 + r2 * r1 + r1**2)

        assert pytest.approx(vol_exp, rel=1e-1) == vol

        fb.add_translation_dof((0,0,1))

        test_matrix = xr.Dataset(coords={
            'rho': rho_water,
            'water_depth': [np.infty],          
            'omega': omega,
            'wave_direction': 0,
            'radiating_dof': list(fb.dofs.keys()),
            })

        solver = BEMSolver()
        data = solver.fill_dataset(test_matrix, 
                                bodies=[fb],
                                hydrostatics=True,
                                mesh=True,
                                wavelength=True,
                                wavenumber=True)
    except ImportError:
        warn('pygmsh and/or meshio is not installed and thus from_meshio has not been tested')

def test_from_meshio_pygmsh_Cylinder():
    try:
        T=0.52
        r1=0.88
        with pygmsh.occ.Geometry() as geom:
                cyl = geom.add_cylinder([0,0,0], [0,0,-T], r1)
                geom.translate(cyl, [0, 0, offset])    
                mesh = geom.generate_mesh(dim=2)
        fb = FloatingBody.from_meshio(mesh)

        vol = fb.mesh.volume
        vol_exp = np.pi*r1**2*T

        assert pytest.approx(vol_exp, rel=1e-1) == vol

        fb.add_translation_dof((0,0,1))

        test_matrix = xr.Dataset(coords={
            'rho': rho_water,
            'water_depth': [np.infty],          
            'omega': omega,
            'wave_direction': 0,
            'radiating_dof': list(fb.dofs.keys()),
            })

        solver = BEMSolver()
        data = solver.fill_dataset(test_matrix, 
                                bodies=[fb],
                                hydrostatics=True,
                                mesh=True,
                                wavelength=True,
                                wavenumber=True)
    except ImportError:
        warn('pygmsh and/or meshio is not installed and thus from_meshio has not been tested')

def test_from_meshio_pygmsh_Sphere():
    try:
        r1=0.88
        with pygmsh.occ.Geometry() as geom:
                sphere = geom.add_ball(center=[0,0,0], radius=r1)
                geom.translate(sphere, [0, 0, offset])    
                mesh = geom.generate_mesh(dim=2)
        fb = FloatingBody.from_meshio(mesh)

        vol = fb.mesh.volume
        vol_exp = 1/2*np.pi*4/3*r1**3

        assert pytest.approx(vol_exp, rel=1e-1) == vol

        fb.add_translation_dof((0,0,1))

        test_matrix = xr.Dataset(coords={
            'rho': rho_water,
            'water_depth': [np.infty],          
            'omega': omega,
            'wave_direction': 0,
            'radiating_dof': list(fb.dofs.keys()),
            })

        solver = BEMSolver()
        data = solver.fill_dataset(test_matrix, 
                                bodies=[fb],
                                hydrostatics=True,
                                mesh=True,
                                wavelength=True,
                                wavenumber=True)
    except ImportError:
        warn('pygmsh and/or meshio is not installed and thus from_meshio has not been tested')


