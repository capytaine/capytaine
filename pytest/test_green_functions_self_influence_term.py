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
"""Tests of the self-influence term of the double layer operators,
computed as the one-sided limit of the integral of the Rankine term wherever a point lies on a face."""

from functools import lru_cache

import pytest
import numpy as np
import capytaine as cpt
import capytaine_test_helpers as helpers


@lru_cache
def flat_square_panel():
    """A single flat square panel, far from the free surface."""
    vertices = np.array([[0.0, 0.0, -2.0], [1.0, 0.0, -2.0], [1.0, 1.0, -2.0], [0.0, 1.0, -2.0]])
    return cpt.Mesh(vertices=vertices, faces=np.array([[0, 1, 2, 3]]))


def full_gradient(points, mesh, **params):
    params.setdefault("free_surface", np.inf)
    _, K = helpers.green_function().evaluate(points, mesh, early_dot_product=False, adjoint_double_layer=True, **params)
    return np.real(K)


def exact_integral_of_rankine_gradient(point, mesh, limit_side):
    """∫∇(1/r) over the first face of mesh, with the given side of the limit (0 for the principal value)."""
    _, int_nabla_G = helpers.green_function().fortran_core.green_rankine.exact_integral_of_rankine(
        point, mesh.vertices[mesh.faces[0], :],
        mesh.faces_centers[0], mesh.faces_normals[0], mesh.faces_areas[0], mesh.faces_radiuses[0],
        True, limit_side)
    return int_nabla_G


@pytest.mark.parametrize("point, expected_fraction", [
    ((0.5, 0.5), 1.0),  # center
    ((0.9, 0.2), 1.0),  # inside, off-center
    ((1.1, 0.5), 0.0),  # outside, closer than the radius
])
def test_in_plane_jump_is_proportional_to_the_fraction_of_the_panel_around_the_point(point, expected_fraction):
    mesh = flat_square_panel()
    point = np.array([*point, -2.0])
    principal_value = exact_integral_of_rankine_gradient(point, mesh, 0.0)
    for limit_side in [-1.0, 1.0]:
        jump = exact_integral_of_rankine_gradient(point, mesh, limit_side) - principal_value
        np.testing.assert_allclose(jump, -2*np.pi*limit_side*expected_fraction*mesh.faces_normals[0], atol=1e-12)


@pytest.mark.parametrize("point, expected_coef", [
    ((0.5, 0.5), 0.5),  # center
    ((0.9, 0.2), 0.5),  # inside, off-center
    ((1.1, 0.5), 0.0),  # outside, closer than the radius
])
def test_self_influence_term_at_points_on_the_panel(point, expected_coef):
    mesh = flat_square_panel()
    point = np.array([*point, -2.0])
    principal_value = -1/(4*np.pi) * exact_integral_of_rankine_gradient(point, mesh, 0.0)
    K = full_gradient(point[None, :], mesh)[:, 0, 0]
    np.testing.assert_allclose(K, principal_value + expected_coef*mesh.faces_normals[0], atol=1e-12)


def test_self_influence_term_is_the_limit_from_the_side_of_the_normal():
    mesh = flat_square_panel()
    eps = 2e-4*mesh.faces_radiuses[0]  # Just above the threshold of the in-plane computation
    for point in [mesh.faces_centers[0], np.array([0.9, 0.2, -2.0])]:
        K_on = full_gradient(point[None, :], mesh)[:, 0, 0]
        K_above = full_gradient((point + eps*mesh.faces_normals[0])[None, :], mesh)[:, 0, 0]
        K_below = full_gradient((point - eps*mesh.faces_normals[0])[None, :], mesh)[:, 0, 0]
        np.testing.assert_allclose(K_on, K_above, atol=1e-3)
        np.testing.assert_allclose(K_above - K_below, mesh.faces_normals[0], atol=1e-3)


def test_potential_finite_at_middle_of_edges():
    mesh = helpers.small_sphere_mesh()
    S, _ = helpers.green_function().evaluate(mesh.middle_edges_waterline, mesh, wavenumber=1.0)
    assert np.all(np.isfinite(S))


@pytest.mark.parametrize("params", [
    dict(free_surface=np.inf),
    dict(wavenumber=1.0),
    dict(wavenumber=1.0, water_depth=10.0),
])
@pytest.mark.parametrize("adjoint_double_layer", [True, False])
@pytest.mark.parametrize("early_dot_product", [True, False])
def test_permuted_collocation_points(params, adjoint_double_layer, early_dot_product):
    mesh = helpers.small_sphere_mesh()
    perm = np.random.default_rng(seed=0).permutation(mesh.nb_faces)
    gf = helpers.green_function()
    kwargs = dict(adjoint_double_layer=adjoint_double_layer, early_dot_product=early_dot_product, **params)
    S_ref, K_ref = gf.evaluate(mesh, mesh, **kwargs)
    S, K = gf.evaluate(mesh.extract_faces(perm), mesh, **kwargs)
    np.testing.assert_allclose(S, S_ref[perm, :], atol=1e-12)
    np.testing.assert_allclose(K, K_ref[..., perm, :], atol=1e-12)


@pytest.mark.parametrize("adjoint_double_layer", [True, False])
def test_hull_rows_of_symmetric_mesh_with_lid(adjoint_double_layer):
    # In the merged mesh including the lid, the faces of the hull and of the lid alternate:
    # the hull faces are not the first faces of the mesh including lid.
    body = helpers.small_sphere_body("rigid", plane_symmetry=True, lid=True)
    mesh_including_lid = body.mesh_including_lid.merged()
    hull = mesh_including_lid.extract_faces(np.where(body.hull_mask)[0])
    gf = helpers.green_function()
    kwargs = dict(wavenumber=1.0, adjoint_double_layer=adjoint_double_layer, early_dot_product=False)
    _, K_ref = gf.evaluate(mesh_including_lid, mesh_including_lid, **kwargs)
    _, K = gf.evaluate(hull, mesh_including_lid, **kwargs)
    np.testing.assert_allclose(K, K_ref[:, body.hull_mask, :], atol=1e-12)
    _, K = gf.evaluate(hull.faces_centers, mesh_including_lid, **kwargs)
    np.testing.assert_allclose(K, K_ref[:, body.hull_mask, :], atol=1e-12)


@pytest.mark.parametrize("wavenumber", [0.0, np.inf])  # low_freq and high_freq singularities, respectively
@pytest.mark.parametrize("adjoint_double_layer", [True, False])
def test_coefficient_one_on_the_free_surface(wavenumber, adjoint_double_layer):
    # Without wave part, the Green function of a flat face on the free surface at its own center is only the
    # Rankine term and its image, whose principal values have no normal component: only the coefficient remains.
    body = helpers.small_sphere_body(lid=True)
    mesh = body.mesh_including_lid
    gf = helpers.green_function()
    _, K = gf.evaluate(mesh, mesh, wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer)
    on_free_surface = np.abs(mesh.faces_centers[:, 2]) < 1e-8
    assert np.any(on_free_surface)
    np.testing.assert_allclose(np.diag(K)[on_free_surface], 1.0, atol=1e-10)


def test_coefficient_one_on_the_sea_bottom():
    depth = 5.0
    mesh = cpt.mesh_parallelepiped(size=(1.0, 1.0, 1.0), center=(0, 0, -depth+0.5), resolution=(2, 2, 2))
    faces = np.where(np.abs(mesh.faces_centers[:, 2] + depth) < 1e-8)[0]
    assert len(faces) > 0
    params = dict(wavenumber=1.0, water_depth=depth, free_surface=0.0)
    _, K = helpers.green_function().evaluate(mesh, mesh, adjoint_double_layer=True, **params)
    diag = np.real(np.diag(K)[faces])
    np.testing.assert_allclose(diag, 1.0, atol=1e-3)
    # Normal gradient at points slightly off these faces, on the side of the normal
    points = mesh.faces_centers[faces] + 2e-4*mesh.faces_radiuses[faces, None]*mesh.faces_normals[faces]
    gradient = full_gradient(points, mesh, **params)
    off_face = np.einsum("ki,ik->i", gradient[:, np.arange(len(faces)), faces], mesh.faces_normals[faces])
    np.testing.assert_allclose(diag, off_face, atol=1e-2)
