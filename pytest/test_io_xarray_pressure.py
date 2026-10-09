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
"""Tests for the pressure-on-hull export added to `fill_dataset`/`assemble_dataset`
via the `keep_details` keyword (`diffraction_pressure`, `Froude_Krylov_pressure`,
`radiation_pressure` dataset variables, and `LinearPotentialFlowResult.pressure_on_hull`).
"""
import pytest

import numpy as np
import xarray as xr

import capytaine as cpt
import capytaine_test_helpers as helpers

from capytaine.meshes.symmetric_meshes import ReflectionSymmetricMesh
from capytaine.meshes.predefined import mesh_horizontal_cylinder, mesh_rectangle


#######################################################################
#                       fill_dataset(keep_details=...)                #
#######################################################################

def test_pressure_variables_present_with_keep_details():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': [1.0, 2.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave'],
    })
    dataset = solver.fill_dataset(test_matrix, sphere, keep_details=True)

    for var in ("diffraction_pressure", "Froude_Krylov_pressure", "radiation_pressure"):
        assert var in dataset
        assert dataset[var].sizes["hull_face"] == sphere.mesh.nb_faces


def test_pressure_variables_absent_by_default():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': [1.0, 2.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave'],
    })
    dataset = solver.fill_dataset(test_matrix, sphere)  # keep_details defaults to False

    for var in ("diffraction_pressure", "Froude_Krylov_pressure", "radiation_pressure"):
        assert var not in dataset


def test_radiation_pressure_has_radiating_dof_dimension():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    test_matrix = xr.Dataset(coords={
        'omega': [1.0], 'radiating_dof': ['Heave', 'Surge'],
    })
    dataset = solver.fill_dataset(test_matrix, sphere, keep_details=True)

    assert "radiating_dof" in dataset["radiation_pressure"].dims
    assert dataset["radiation_pressure"].sizes["radiating_dof"] == 2

    heave_pressure = dataset["radiation_pressure"].sel(radiating_dof="Heave").values
    surge_pressure = dataset["radiation_pressure"].sel(radiating_dof="Surge").values
    assert not np.allclose(heave_pressure, surge_pressure)


def test_failed_results_do_not_crash_pressure_export():
    sphere = helpers.small_sphere_body(["Heave"])
    broken_bem_solver = helpers.broken_bem_solver()
    test_matrix = xr.Dataset(coords={
        "wavenumber": np.linspace(0.1, 5.0, 5), "wave_direction": [0.0], "radiating_dof": ["Heave"],
    })
    dataset = broken_bem_solver.fill_dataset(test_matrix, sphere, keep_details=True)

    assert np.any(np.isnan(dataset.added_mass.values))
    assert np.any(np.isnan(dataset["radiation_pressure"].values))
    # Froude-Krylov pressure is computed analytically from the problem alone,
    # so it should stay finite even for the wavenumbers where the solve failed.
    assert np.all(np.isfinite(dataset["Froude_Krylov_pressure"].values))


def test_pressure_ignored_for_infinite_free_surface():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb = cpt.RadiationProblem(body=sphere, free_surface=np.inf, radiating_dof="Heave")
    res = solver.solve(pb, keep_details=True)
    dataset = cpt.assemble_dataset([res])
    assert "radiation_pressure" not in dataset


def test_keep_details_forced_by_kochin_export():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    test_matrix = xr.Dataset(coords={
        'omega': [1.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave'],
        'theta': np.linspace(0, 2*np.pi, 5),
    })
    dataset = solver.fill_dataset(test_matrix, sphere, keep_details=False)
    assert "radiation_pressure" in dataset
    assert "diffraction_pressure" in dataset


#######################################################################
#                  LinearPotentialFlowResult.pressure_on_hull          #
#######################################################################

def test_pressure_on_hull_matches_hull_mask():
    sphere = helpers.small_sphere_body(["Heave"])
    solver = helpers.solver()
    pb = cpt.RadiationProblem(body=sphere, omega=1.0, radiating_dof="Heave")
    res = solver.solve(pb, keep_details=True)
    np.testing.assert_array_equal(res.pressure_on_hull, res.pressure[sphere.hull_mask])
    assert res.pressure_on_hull.shape == (sphere.mesh.nb_faces,)


def test_pressure_on_hull_with_symmetric_mesh_and_lid():
    mesh = mesh_horizontal_cylinder(reflection_symmetry=True).immersed_part()
    lid_mesh = ReflectionSymmetricMesh(
            mesh_rectangle(size=(1.0, 10.0), faces_max_radius=0.5, center=(0, 0.5, -0.05)),
            plane=mesh.plane
            )
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = helpers.solver()
    res = solver.solve(pb, keep_details=True)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb, keep_details=True)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])
    np.testing.assert_allclose(
            np.sort_complex(res.pressure_on_hull),
            np.sort_complex(ref_res.pressure_on_hull),
            rtol=1e-3,
            )
