#!/usr/bin/env python
# coding: utf-8

import pytest
import numpy as np
import xarray as xr
import capytaine as cpt


builtin_methods = ['GL2']
quadpy_methods = ['sommariva_05']

@pytest.mark.parametrize('method', builtin_methods + quadpy_methods)
def test_quadrature_points_are_coplanar(method):
    # Check that all quadrature points are within the plane containing the panel
    mesh = cpt.mesh_sphere()
    if method not in builtin_methods:
        quadpy = pytest.importorskip("quadpy", reason="Quadpy not installed, test skipped.")
        method = quadpy.c2.schemes[method]()
    mesh.compute_quadrature(method)
    for i_face in range(mesh.nb_faces):
        A, B, C, D = mesh.vertices[mesh.faces[i_face, :], :]
        for j_quad_point in range(mesh.quadrature_points[0].shape[1]):
            X = mesh.quadrature_points[0][i_face, j_quad_point, :]
            assert np.isclose(np.dot(np.cross(B-A, C-A), X-A), 0.0, atol=1e-8)


@pytest.mark.parametrize('method', builtin_methods + quadpy_methods)
def test_area(method):
    # Check that quadrature weights sum to the area of the panel
    mesh = cpt.mesh_sphere()
    if method not in builtin_methods:
        quadpy = pytest.importorskip("quadpy", reason="Quadpy not installed, test skipped.")
        method = quadpy.c2.schemes[method]()
    mesh.compute_quadrature(method)
    for i_face in range(mesh.nb_faces):
        assert np.isclose(np.sum(mesh.quadrature_points[1][i_face, :]), mesh.faces_areas[i_face], rtol=1e-2)


@pytest.mark.parametrize('method', builtin_methods + quadpy_methods)
def test_resolution(method):
    mesh = cpt.mesh_horizontal_cylinder(
        length=5.0, radius=1.0,
        center=(0, 0, 0),
        resolution=(2, 5, 10),
    ).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Heave")

    test_matrix = xr.Dataset(coords={
        "omega": np.linspace(0.5, 3.0, 2),
        "radiating_dof": ["Heave"],
    })

    solver = cpt.BEMSolver()

    data_0 = solver.fill_dataset(test_matrix, [body], mesh=True)

    if method not in builtin_methods:
        quadpy = pytest.importorskip("quadpy", reason="Quadpy not installed, test skipped.")
        body.mesh.compute_quadrature(quadpy.c2.schemes[method]())
    else:
        body.mesh.compute_quadrature(method)
    data_1 = solver.fill_dataset(test_matrix, [body], mesh=True)

    assert np.allclose(data_0["added_mass"].data, data_1["added_mass"].data, rtol=1e-1)
