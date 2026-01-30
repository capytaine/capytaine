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
    assert "kochin_radiation" not in kds_diff
    assert kds_diff.sizes == ({'wavelength': 1, 'wave_direction': 1, 'theta': 3})

def test_kochin_array_radiation():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_rad = solver.solve(pb_rad)
    kds_rad = kochin_data_array([res_rad], np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin_radiation" in kds_rad
    assert kds_rad.sizes == ({'wavelength': 1, 'radiating_dof': 1, 'theta': 3})

def test_kochin_array_both_radiation_and_diffraction():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0)
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_both = solver.solve_all([pb_rad, pb_diff])
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin_radiation" in kds_both
    assert kds_both.sizes == ({'wavelength': 1, 'wave_direction': 1, 'radiating_dof': 1, 'theta': 3})

def test_kochin_array_diffraction_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0, forward_speed=u) for u in [0.0, 1.0]]
    res_diff = solver.solve_all(pb_diff)
    kds_diff = kochin_data_array(res_diff, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_diff
    assert "kochin_radiation" not in kds_diff
    assert kds_diff.sizes == ({'wavelength': 1, 'wave_direction': 1, 'forward_speed': 2, 'theta': 3})

def test_kochin_array_radiation_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", forward_speed=u, wave_direction=pi) for u in [0.0, 1.0]]
    res_rad = solver.solve_all(pb_rad)
    kds_rad = kochin_data_array(res_rad, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin_radiation" in kds_rad
    assert kds_rad.sizes == ({'wavelength': 1, 'radiating_dof': 1, 'forward_speed': 2, 'theta': 3})

def test_kochin_array_both_radiation_and_diffraction_with_forward_speed():
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    res_both = solver.solve_all(pb_rad + pb_diff)
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin_radiation" in kds_both
    assert kds_both.sizes == ({'wavelength': 1, 'wave_direction': 1, 'radiating_dof': 1, 'forward_speed': 2, 'theta': 3})

def test_fill_dataset_with_kochin_functions():
    mesh = cpt.mesh_sphere().immersed_part()
    sphere = cpt.FloatingBody(mesh=mesh)
    sphere.add_translation_dof(name="Heave")
    solver = cpt.BEMSolver()

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': [0.0, pi/2],
        'radiating_dof': ["Heave"],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin_radiation' in ds
    assert 'kochin_diffraction' not in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin_radiation'].sel(radiating_dof="Heave", theta=0.0),
                      ds['kochin_radiation'].sel(radiating_dof="Heave", theta=pi/2))

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'radiating_dof': ["Heave"],
        'wave_direction': [-pi/2, 0.0],
        'theta': [0.0, pi/2],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin_radiation' in ds
    assert 'kochin_diffraction' in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin_diffraction'].sel(wave_direction=-pi/2, theta=0.0),
                      ds['kochin_diffraction'].sel(wave_direction=0.0, theta=pi/2))
    
def test_kochin_mesh_lid():
    mesh = cpt.mesh_sphere().immersed_part()
    lid_mesh = mesh.generate_lid()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), lid_mesh=lid_mesh)
    problem = cpt.RadiationProblem(body=body, radiating_dof='Surge', omega=0.3)
    body_without_lid = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    body_with_lid = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), lid_mesh=lid_mesh)
    omega = 0.3
    problem_without_lid = cpt.RadiationProblem(body=body_without_lid, radiating_dof='Surge', omega=omega)
    problem_with_lid = cpt.RadiationProblem(body=body_with_lid, radiating_dof='Surge', omega=omega)
    solver = cpt.BEMSolver()
    res = solver.solve(problem)
    res_without_lid = solver.solve(problem_without_lid)
    res_with_lid = solver.solve(problem_with_lid)
    theta = 0.2
    H = cpt.post_pro.compute_kochin(res, theta)
    H_without_lid = cpt.post_pro.compute_kochin(res_without_lid, theta)
    H_with_lid = cpt.post_pro.compute_kochin(res_with_lid, theta)
    assert np.isclose(H_without_lid, H_with_lid, atol=1e-05)
