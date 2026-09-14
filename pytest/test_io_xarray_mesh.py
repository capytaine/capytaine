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
"""Tests for the mesh/dofs export added to `fill_dataset`/`assemble_dataset`
via the `mesh=True` keyword (`mesh_vertices`, `mesh_faces_center`,
`lid_mesh_vertices`, `lid_mesh_faces_center`, `dof_motions` and
`dof_gradient_of_motions` dataset variables).
"""
import pytest

import numpy as np
import xarray as xr

import capytaine as cpt
from capytaine.meshes.predefined import mesh_horizontal_cylinder


@pytest.fixture
def sphere():
    sphere = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, -2), radius=1.0, resolution=(4, 4)),
            dofs=cpt.rigid_body_dofs(),
            name="sphere",
            )
    return sphere


@pytest.fixture
def solver():
    return cpt.BEMSolver()


def test_mesh_vertices_and_faces_center(sphere, solver):
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, sphere, mesh=True, hydrostatics=False)

    nb_faces = sphere.mesh.nb_faces
    assert dataset.sizes["hull_face"] == nb_faces
    assert dataset["mesh_vertices"].dims == ("hull_face", "vertices_of_face", "space_coordinate")
    assert dataset["mesh_vertices"].shape[0] == nb_faces
    assert dataset["mesh_faces_center"].dims == ("hull_face", "space_coordinate")
    np.testing.assert_allclose(dataset["mesh_faces_center"].values, sphere.mesh.faces_centers)
    assert "quadrature_method" in dataset.coords
    assert "lid_mesh_vertices" not in dataset
    assert "lid_mesh_faces_center" not in dataset


def test_mesh_export_shares_hull_face_dimension_with_pressure(sphere, solver):
    # mesh_faces_center should let an external reader recover the position of
    # each `hull_face` index used by the pressure variables.
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, sphere, mesh=True, keep_details=True, hydrostatics=False)
    assert dataset["mesh_faces_center"].sizes["hull_face"] == dataset["radiation_pressure"].sizes["hull_face"]


def test_dof_motions_matches_dof_evaluation(sphere, solver):
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave', 'Surge']})
    dataset = solver.fill_dataset(test_matrix, sphere, mesh=True, hydrostatics=False)

    assert "influenced_dof" in dataset["dof_motions"].dims
    np.testing.assert_allclose(
            dataset["dof_motions"].sel(influenced_dof="Heave").values,
            sphere.dofs["Heave"].evaluate_motion(sphere.mesh),
            )
    np.testing.assert_allclose(
            dataset["dof_motions"].sel(influenced_dof="Surge").values,
            sphere.dofs["Surge"].evaluate_motion(sphere.mesh),
            )


def test_dof_gradient_of_motions(sphere, solver):
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave', 'Roll']})
    dataset = solver.fill_dataset(test_matrix, sphere, mesh=True, hydrostatics=False)

    grad = dataset["dof_gradient_of_motions"]
    assert grad.dims == ("influenced_dof", "hull_face", "space_coordinate", "gradient")
    assert not np.issubdtype(grad.dtype, np.complexfloating)

    # Heave is a rigid translation: the gradient of its motion is exactly zero everywhere.
    np.testing.assert_allclose(grad.sel(influenced_dof="Heave").values, 0.0)

    # Roll is a rotation: the gradient should match evaluate_gradient_of_motion directly.
    np.testing.assert_allclose(
            grad.sel(influenced_dof="Roll").values,
            sphere.dofs["Roll"].evaluate_gradient_of_motion(sphere.mesh),
            )


def test_mesh_export_with_raw_array_dof(sphere, solver):
    sphere.dofs = {"Heave": sphere.dofs["Heave"], "LegacyCustom": np.ones((sphere.mesh.nb_faces, 3))}
    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave', 'LegacyCustom']})
    dataset = solver.fill_dataset(test_matrix, sphere, mesh=True, hydrostatics=False)

    np.testing.assert_allclose(dataset["dof_motions"].sel(influenced_dof="LegacyCustom").values, 1.0)
    # No analytic gradient is available for a plain-array dof: filled with NaN,
    # and it should not force the whole variable to a complex dtype.
    custom_grad = dataset["dof_gradient_of_motions"].sel(influenced_dof="LegacyCustom")
    assert not np.issubdtype(custom_grad.dtype, np.complexfloating)
    assert np.all(np.isnan(custom_grad.values))
    heave_grad = dataset["dof_gradient_of_motions"].sel(influenced_dof="Heave")
    assert np.all(np.isfinite(heave_grad.values))


def test_mesh_export_with_lid(solver):
    mesh = cpt.mesh_parallelepiped(center=(0, 0, -1.0), size=(2.0, 2.0, 2.0))
    hull_mesh, lid_mesh = mesh.extract_lid()
    body = cpt.FloatingBody(mesh=hull_mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs(only=["Heave"]))

    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, body, mesh=True, hydrostatics=False)

    assert dataset.sizes["lid_face"] == lid_mesh.nb_faces
    assert dataset["lid_mesh_vertices"].dims == ("lid_face", "vertices_of_face", "space_coordinate")
    np.testing.assert_allclose(dataset["lid_mesh_faces_center"].values, lid_mesh.faces_centers)
    # The hull-only dimension should not be affected by the lid.
    assert dataset.sizes["hull_face"] == body.mesh.nb_faces


def test_mesh_export_with_multibody_dof_on_submesh(solver):
    body1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, -2), radius=1.0, resolution=(4, 4)),
            dofs=cpt.rigid_body_dofs(only=["Heave"]),
            name="body1",
            )
    body2 = body1.translated_x(10.0, name="body2")
    both = body1 + body2

    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': list(both.dofs.keys())})
    dataset = solver.fill_dataset(test_matrix, [both], mesh=True, hydrostatics=False)

    assert dataset.sizes["hull_face"] == both.mesh.nb_faces
    for dof_name, dof in both.dofs.items():
        np.testing.assert_allclose(
                dataset["dof_motions"].sel(influenced_dof=dof_name).values,
                dof.evaluate_motion(both.mesh),
                )


def test_mesh_export_with_symmetric_mesh(solver):
    sym_mesh = mesh_horizontal_cylinder(reflection_symmetry=True).immersed_part()
    sym_body = cpt.FloatingBody(mesh=sym_mesh, dofs=cpt.rigid_body_dofs(only=["Heave"]), name="sym_body")

    test_matrix = xr.Dataset(coords={'omega': [1.0], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, sym_body, mesh=True, hydrostatics=False)

    assert dataset.sizes["hull_face"] == sym_mesh.nb_faces
    assert dataset.coords["quadrature_method"].values.item() == str(sym_mesh.quadrature_method)
    np.testing.assert_allclose(dataset["mesh_faces_center"].values, sym_mesh.faces_centers)
    np.testing.assert_allclose(
            dataset["dof_motions"].sel(influenced_dof="Heave").values,
            sym_body.dofs["Heave"].evaluate_motion(sym_body.mesh),
            )

    # The symmetry itself is not preserved: the exported mesh is numerically
    # identical to that of the equivalent non-symmetric (merged) body, and
    # the dataset carries no indication that the original mesh was symmetric.
    merged_body = cpt.FloatingBody(mesh=sym_mesh.merged(), dofs=cpt.rigid_body_dofs(only=["Heave"]), name="merged_body")
    merged_dataset = solver.fill_dataset(test_matrix, merged_body, mesh=True, hydrostatics=False)
    np.testing.assert_allclose(dataset["mesh_vertices"].values, merged_dataset["mesh_vertices"].values)
    np.testing.assert_allclose(dataset["mesh_faces_center"].values, merged_dataset["mesh_faces_center"].values)
