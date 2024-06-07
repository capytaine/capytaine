"""Tests for the mesh import of mesh from external libraries"""

import os
import pytest

import numpy as np
import xarray as xr

from capytaine.io.mesh_writers import write_STL
from capytaine.io.mesh_loaders import load_STL
import capytaine as cpt

#################################################################

def generate_pygmsh_wavebot():
    pygmsh = pytest.importorskip("pygmsh", reason="PyGMSH not installed, test skipped")
    offset = 1e-2 # small offset so you can clip at z=0
    T1 = 0.16
    T2 = 0.37
    r1 = 0.88
    r2 = 0.35
    with pygmsh.occ.Geometry() as geom:
        cyl = geom.add_cylinder([0,0,0], [0,0,-T1], r1)
        cone = geom.add_cone([0,0,-T1], [0,0,-T2], r1, r2)
        geom.translate(cyl, [0, 0, offset])
        geom.translate(cone, [0, 0, offset])
        geom.boolean_union([cyl, cone])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = np.pi*T1*r1**2 + 1/3 * np.pi * T2 * (r2**2 + r2 * r1 + r1**2)
    return (mesh, vol_exp)

def generate_pygmsh_cylinder():
    pygmsh = pytest.importorskip("pygmsh", reason="PyGMSH not installed, test skipped")
    offset = 1e-2 # small offset so you can clip at z=0
    T = 0.52
    r1 = 0.88
    with pygmsh.occ.Geometry() as geom:
        cyl = geom.add_cylinder([0,0,0], [0,0,-T], r1)
        geom.translate(cyl, [0, 0, offset])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = np.pi*r1**2*T
    return (mesh, vol_exp)

def generate_pygmsh_sphere():
    pygmsh = pytest.importorskip("pygmsh", reason="PyGMSH not installed, test skipped")
    offset = 1e-2 # small offset so you can clip at z=0
    r1 = 0.88
    with pygmsh.occ.Geometry() as geom:
        sphere = geom.add_ball(center=[0,0,0], radius=r1)
        geom.translate(sphere, [0, 0, offset])
        mesh = geom.generate_mesh(dim=2)
    vol_exp = 1/2*np.pi*4/3*r1**3
    return (mesh, vol_exp)


@pytest.mark.parametrize("generate_pygmsh",
                         [generate_pygmsh_cylinder, generate_pygmsh_sphere, generate_pygmsh_wavebot])
def test_from_meshio_pygmsh(generate_pygmsh):
    pytest.importorskip("pygmsh", reason="PyGMSH not installed, test skipped")
    mesh, vol_exp = generate_pygmsh()
    mesh = cpt.load_mesh(mesh)
    assert mesh.immersed_part().volume == pytest.approx(vol_exp, rel=1e-1)

    fb = cpt.FloatingBody(mesh)
    fb.add_translation_dof(name="Heave")
    test_matrix = xr.Dataset(coords={
        'omega': [np.pi],
        'wave_direction': [0],
        'radiating_dof': list(fb.dofs.keys()),
        })
    solver = cpt.BEMSolver()
    solver.fill_dataset(test_matrix, bodies=fb)

#################################################################

def test_STL(tmp_path):
    pytest.importorskip("vtk", reason="VTK not installed, test skipped.")
    mesh = cpt.mesh_sphere()
    filepath = tmp_path / "test.stl"
    write_STL(filepath, mesh.vertices, mesh.faces)
    reloaded_mesh = load_STL(str(filepath), name="Bla")

    assert reloaded_mesh.name == "Bla"
    assert np.allclose(mesh.vertices, reloaded_mesh.vertices)


def test_STL_with_uppercase_file_extension(tmp_path):
    pytest.importorskip("vtk", reason="VTK not installed, test skipped.")
    mesh = cpt.mesh_sphere()
    filepath = tmp_path / "test.STL"
    write_STL(filepath, mesh.vertices, mesh.faces)
    reloaded_mesh = cpt.load_mesh(str(filepath), name="Bla")
    # Should detect that the file is stl using the file extension

    assert reloaded_mesh.name == "Bla"
    assert np.allclose(mesh.vertices, reloaded_mesh.vertices)


def test_STL_with_meshio(tmp_path):
    pytest.importorskip("vtk", reason="VTK not installed, test skipped.")
    pytest.importorskip("meshio", reason="Meshio not installed, test skipped")
    from capytaine.io.meshio import load_from_meshio
    mesh, _ = generate_pygmsh_cylinder()
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

#################################################################

def test_MED_file():
    pytest.importorskip("h5py", reason="h5py not installed, test skipped")
    mesh = cpt.load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/barge.med"))
    assert mesh.nb_faces == 187

def test_MSH2_file():
    mesh = cpt.load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder2.msh"))
    assert mesh.nb_faces == 64

def test_MSH4_file():
    pytest.importorskip("meshio", reason="meshio not installed, test skipped")
    mesh = cpt.load_mesh(os.path.join(os.path.dirname(__file__), "mesh_files_examples/cylinder4.msh"))
    assert mesh.nb_faces == 64
