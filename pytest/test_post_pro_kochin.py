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
"""Test related to the computation of Kochin functions."""

import numpy as np
from numpy import pi
import pandas as pd
import xarray as xr
import pytest

import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.io.xarray import kochin_data_array


def test_kochin_array_diffraction():
    body = helpers.small_sphere_body()
    solver = helpers.solver()
    pb_diff = cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0)
    res_diff = solver.solve(pb_diff)
    kds_diff = kochin_data_array([res_diff], np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_diff
    assert "kochin_radiation" not in kds_diff
    assert kds_diff.sizes == ({'wavelength': 1, 'wave_direction': 1, 'theta': 3})

def test_kochin_array_radiation():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_rad = solver.solve(pb_rad)
    kds_rad = kochin_data_array([res_rad], np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin_radiation" in kds_rad
    assert kds_rad.sizes == ({'wavelength': 1, 'radiating_dof': 1, 'theta': 3})

def test_kochin_array_both_radiation_and_diffraction():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb_diff = cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0)
    pb_rad = cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave")
    res_both = solver.solve_all([pb_rad, pb_diff])
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin_radiation" in kds_both
    assert kds_both.sizes == ({'wavelength': 1, 'wave_direction': 1, 'radiating_dof': 1, 'theta': 3})

def test_kochin_array_diffraction_with_forward_speed():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=0.0, forward_speed=u) for u in [0.0, 1.0]]
    res_diff = solver.solve_all(pb_diff)
    kds_diff = kochin_data_array(res_diff, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_diff
    assert "kochin_radiation" not in kds_diff
    assert kds_diff.sizes == ({'wavelength': 1, 'wave_direction': 1, 'forward_speed': 2, 'theta': 3})

def test_kochin_array_radiation_with_forward_speed():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", forward_speed=u, wave_direction=pi) for u in [0.0, 1.0]]
    res_rad = solver.solve_all(pb_rad)
    kds_rad = kochin_data_array(res_rad, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" not in kds_rad
    assert "kochin_radiation" in kds_rad
    assert kds_rad.sizes == ({'wavelength': 1, 'radiating_dof': 1, 'forward_speed': 2, 'theta': 3})

def test_kochin_array_both_radiation_and_diffraction_with_forward_speed():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb_diff = [cpt.DiffractionProblem(body=body, wavelength=2.0, wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    pb_rad = [cpt.RadiationProblem(body=body, wavelength=2.0, radiating_dof="Heave", wave_direction=pi, forward_speed=u) for u in [0.0, 1.0]]
    res_both = solver.solve_all(pb_rad + pb_diff)
    kds_both = kochin_data_array(res_both, np.linspace(0.0, np.pi, 3))
    assert "kochin_diffraction" in kds_both
    assert "kochin_radiation" in kds_both
    assert kds_both.sizes == ({'wavelength': 1, 'wave_direction': 1, 'radiating_dof': 1, 'forward_speed': 2, 'theta': 3})

def test_fill_dataset_with_kochin_functions():
    body = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': [0.0, pi/2],
        'radiating_dof': ["Heave"],
    })
    ds = solver.fill_dataset(test_matrix, body)
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
    ds = solver.fill_dataset(test_matrix, body)
    assert 'theta' in ds.coords
    assert 'kochin_radiation' in ds
    assert 'kochin_diffraction' in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin_diffraction'].sel(wave_direction=-pi/2, theta=0.0),
                      ds['kochin_diffraction'].sel(wave_direction=0.0, theta=pi/2))

def test_kochin_mesh_lid():
    body_without_lid = helpers.small_sphere_body(lid=False)
    body_with_lid = helpers.small_sphere_body(lid=True)
    omega = 0.3
    problem_without_lid = cpt.RadiationProblem(body=body_without_lid, radiating_dof='Surge', omega=omega)
    problem_with_lid = cpt.RadiationProblem(body=body_with_lid, radiating_dof='Surge', omega=omega)
    solver = helpers.solver()
    res_without_lid = solver.solve(problem_without_lid)
    res_with_lid = solver.solve(problem_with_lid)
    theta = 0.2
    H_without_lid = cpt.post_pro.compute_kochin(res_without_lid, theta)
    H_with_lid = cpt.post_pro.compute_kochin(res_with_lid, theta)
    assert np.isclose(H_without_lid, H_with_lid, atol=1e-05)
