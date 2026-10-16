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
"""Tests for the sources (indirect method) and potential (direct method) stored in the dataset
by `fill_dataset(..., keep_details=True)`."""
import numpy as np
import xarray as xr

import capytaine as cpt
import capytaine_test_helpers as helpers

from capytaine.io.xarray import separate_complex_values, merge_complex_values, save_dataset_as_netcdf

SOURCES = ["diffraction_sources", "radiation_sources", "diffraction_lid_sources", "radiation_lid_sources"]
POTENTIALS = ["diffraction_potential", "radiation_potential", "diffraction_lid_potential", "radiation_lid_potential"]
POINTS = np.array([[0.3, 0.2, -1.5], [2.0, 1.0, -0.5], [-1.0, 3.0, -2.0]])


def _test_matrix():
    return xr.Dataset(coords={
        'omega': [1.0, 2.0], 'wave_direction': [0.0, np.pi/3], 'radiating_dof': ['Heave', 'Surge'],
    })


def _body(lid):
    return helpers.small_sphere_body(["Heave", "Surge"], lid=lid)


def _all_face_values(body, hull, lid):
    """Put back hull and lid values in the order of `mesh_including_lid`."""
    values = np.zeros(body.mesh_including_lid.nb_faces, dtype=complex)
    values[body.hull_mask] = hull
    if lid is not None:
        values[~body.hull_mask] = lid
    return values


def _check_roundtrip_of_velocity(dataset, body, solver, lid):
    for omega in dataset.omega.values:
        for beta in dataset.wave_direction.values:
            sources = _all_face_values(
                body,
                dataset.diffraction_sources.sel(omega=omega, wave_direction=beta).values,
                dataset.diffraction_lid_sources.sel(omega=omega, wave_direction=beta).values if lid else None,
            )
            problem = cpt.DiffractionProblem(body=body, omega=omega, wave_direction=beta)
            ref = solver.solve(problem, keep_details=True)
            from_dataset = problem.make_results_container(sources=sources)
            assert np.allclose(solver.compute_velocity(POINTS, from_dataset), solver.compute_velocity(POINTS, ref))
        for dof in dataset.radiating_dof.values:
            sources = _all_face_values(
                body,
                dataset.radiation_sources.sel(omega=omega, radiating_dof=dof).values,
                dataset.radiation_lid_sources.sel(omega=omega, radiating_dof=dof).values if lid else None,
            )
            problem = cpt.RadiationProblem(body=body, omega=omega, radiating_dof=dof)
            ref = solver.solve(problem, keep_details=True)
            from_dataset = problem.make_results_container(sources=sources)
            assert np.allclose(solver.compute_velocity(POINTS, from_dataset), solver.compute_velocity(POINTS, ref))


def test_sources_sizes_without_lid():
    body = _body(lid=False)
    dataset = helpers.solver().fill_dataset(_test_matrix(), body, keep_details=True, hydrostatics=False)
    assert dataset.diffraction_sources.dims == ("omega", "wave_direction", "hull_face")
    assert dataset.radiation_sources.dims == ("omega", "radiating_dof", "hull_face")
    assert dataset.diffraction_sources.sizes["hull_face"] == body.mesh.nb_faces
    assert dataset.radiation_sources.sizes["hull_face"] == body.mesh.nb_faces
    assert "diffraction_lid_sources" not in dataset
    assert "radiation_lid_sources" not in dataset
    for var in POTENTIALS:
        assert var not in dataset


def test_sources_sizes_with_lid():
    body = _body(lid=True)
    dataset = helpers.solver().fill_dataset(_test_matrix(), body, keep_details=True, hydrostatics=False)
    assert dataset.diffraction_sources.sizes["hull_face"] == body.mesh.nb_faces
    assert dataset.radiation_sources.sizes["hull_face"] == body.mesh.nb_faces
    assert dataset.diffraction_lid_sources.dims == ("omega", "wave_direction", "lid_face")
    assert dataset.radiation_lid_sources.dims == ("omega", "radiating_dof", "lid_face")
    assert dataset.diffraction_lid_sources.sizes["lid_face"] == body.lid_mesh.nb_faces
    assert dataset.radiation_lid_sources.sizes["lid_face"] == body.lid_mesh.nb_faces


def test_sources_absent_by_default():
    dataset = helpers.solver().fill_dataset(_test_matrix(), _body(lid=True), hydrostatics=False)
    for var in SOURCES + POTENTIALS:
        assert var not in dataset


def test_sources_squeezed_optional_dims():
    body = _body(lid=False)
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave']})
    dataset = helpers.solver().fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    assert set(dataset.diffraction_sources.dims) == {"omega", "wave_direction", "hull_face"}
    assert set(dataset.radiation_sources.dims) == {"omega", "radiating_dof", "hull_face"}


def test_sources_roundtrip_without_lid():
    body = _body(lid=False)
    solver = helpers.solver()
    dataset = solver.fill_dataset(_test_matrix(), body, keep_details=True, hydrostatics=False)
    _check_roundtrip_of_velocity(dataset, body, solver, lid=False)


def test_sources_roundtrip_with_lid():
    body = _body(lid=True)
    solver = helpers.solver()
    dataset = solver.fill_dataset(_test_matrix(), body, keep_details=True, hydrostatics=False)
    _check_roundtrip_of_velocity(dataset, body, solver, lid=True)


def test_sources_roundtrip_after_separate_complex_values_and_netcdf(tmp_path):
    body = _body(lid=True)
    solver = helpers.solver()
    dataset = solver.fill_dataset(_test_matrix(), body, keep_details=True, hydrostatics=False)

    separated = separate_complex_values(dataset)
    for var in SOURCES:
        assert separated[var].dims[0] == "complex" and not np.iscomplexobj(separated[var])
    assert np.allclose(merge_complex_values(separated).diffraction_lid_sources, dataset.diffraction_lid_sources)

    path = tmp_path / "dataset.nc"
    save_dataset_as_netcdf(path, separated)
    loaded = merge_complex_values(xr.load_dataset(path))
    for var in SOURCES:
        assert loaded[var].dims == dataset[var].dims
        assert np.allclose(loaded[var].values, dataset[var].values)
    _check_roundtrip_of_velocity(loaded, body, solver, lid=True)


def test_direct_method_stores_potential_not_sources():
    body = _body(lid=True)
    test_matrix = xr.Dataset(coords={'omega': [1.0, 2.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave', 'Surge']})
    solver = helpers.solver(method="direct")
    dataset = solver.fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    for var in SOURCES:
        assert var not in dataset
    assert dataset.diffraction_potential.dims == ("omega", "wave_direction", "hull_face")
    assert dataset.radiation_potential.dims == ("omega", "radiating_dof", "hull_face")
    assert dataset.diffraction_lid_potential.dims == ("omega", "wave_direction", "lid_face")
    assert dataset.radiation_lid_potential.dims == ("omega", "radiating_dof", "lid_face")
    assert dataset.radiation_potential.sizes["hull_face"] == body.mesh.nb_faces
    assert dataset.radiation_lid_potential.sizes["lid_face"] == body.lid_mesh.nb_faces

    for omega in test_matrix.omega.values:
        res = solver.solve(cpt.RadiationProblem(body=body, omega=omega, radiating_dof="Surge"), keep_details=True)
        assert np.allclose(dataset.radiation_potential.sel(omega=omega, radiating_dof="Surge").values,
                           res.potential[body.hull_mask])
        assert np.allclose(dataset.radiation_lid_potential.sel(omega=omega, radiating_dof="Surge").values,
                           res.potential[~body.hull_mask])


def test_direct_method_potential_without_lid():
    body = _body(lid=False)
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'wave_direction': [0.0], 'radiating_dof': ['Heave']})
    dataset = helpers.solver(method="direct").fill_dataset(test_matrix, body, keep_details=True, hydrostatics=False)
    assert "diffraction_potential" in dataset and "radiation_potential" in dataset
    assert "diffraction_lid_potential" not in dataset
    assert "radiation_lid_potential" not in dataset


def test_variables_attributes():
    dataset = helpers.solver().fill_dataset(_test_matrix(), _body(lid=True), keep_details=True, hydrostatics=False)
    for var in SOURCES:
        assert "long_name" in dataset[var].attrs
