"""Test related to the computation of Kochin functions."""

import numpy as np
from numpy import pi
import pandas as pd
import xarray as xr
import pytest

import capytaine as cpt
from capytaine.io.xarray import kochin_data_array


def test_kochin_array_diffraction():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0)
    res_diff = solver.solve(pb_diff)
    kds_diff = kochin_data_array([res_diff], np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_diff
    assert "kochin" not in kds_diff
    assert len(kds_diff.dims) == 3

def test_kochin_array_radiation():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_rad = solver.solve(pb_rad)
    kds_rad = kochin_data_array([res_rad], np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin" in kds_rad
    assert len(kds_rad.dims) == 3

def test_kochin_array_both_radiation_and_diffraction():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0)
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_both = solver.solve_all([pb_rad, pb_diff])
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin" in kds_both
    assert len(kds_both.dims) == 4

def test_kochin_array_diffraction_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0, forward_speed=u) for u in [0.0, 1.0]]
    res_diff = solver.solve_all(pb_diff)
    kds_diff = kochin_data_array(res_diff, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_diff
    assert "kochin" not in kds_diff
    assert len(kds_diff.dims) == 4

def test_kochin_array_radiation_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", forward_speed=u, wave_direction=pi) for u in [0.0, 1.0]]
    res_rad = solver.solve_all(pb_rad)
    kds_rad = kochin_data_array(res_rad, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin" in kds_rad
    assert len(kds_rad.dims) == 4

def test_kochin_array_both_radiation_and_diffraction_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    res_both = solver.solve_all(pb_rad + pb_diff)
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin" in kds_both
    assert len(kds_both.dims) == 5

def test_fill_dataset_with_kochin_functions():
    sphere = cpt.Sphere(clip_free_surface=True)
    sphere.add_translation_dof(name="Heave")
    solver = cpt.BEMSolver()

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': [0.0, pi/2],
        'radiating_dof': ["Heave"],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin' in ds
    assert 'kochin_diffraction' not in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin'].sel(radiating_dof="Heave", theta=0.0),
                      ds['kochin'].sel(radiating_dof="Heave", theta=pi/2))

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'radiating_dof': ["Heave"],
        'wave_direction': [-pi/2, 0.0],
        'theta': [0.0, pi/2],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin' in ds
    assert 'kochin_diffraction' in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin_diffraction'].sel(wave_direction=-pi/2, theta=0.0),
                      ds['kochin_diffraction'].sel(wave_direction=0.0, theta=pi/2))
