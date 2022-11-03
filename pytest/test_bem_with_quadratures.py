#!/usr/bin/env python
# coding: utf-8

import pytest
import numpy as np
import xarray as xr
import capytaine as cpt

try:
    import quadpy
except ImportError:
    quadpy = None


@pytest.mark.skipif(quadpy is None, reason="quadpy is not installed")
def test_quadrature_points_are_coplanar():
    # Check that all quadrature points are within the plane containing the panel
    mesh = cpt.Sphere().mesh
    for method in [quadpy.c2.schemes['sommariva_05'](), quadpy.c2.schemes["stroud_c2_7_4"]()]:
        mesh.compute_quadrature(method)
        for i_face in range(mesh.nb_faces):
            A, B, C, D = mesh.vertices[mesh.faces[i_face, :], :]
            for j_quad_point in range(mesh.quadrature_points[0].shape[1]):
                X = mesh.quadrature_points[0][i_face, j_quad_point, :]
                assert np.isclose(np.dot(np.cross(B-A, C-A), X-A), 0.0, atol=1e-8)


@pytest.mark.skipif(quadpy is None, reason="quadpy is not installed")
def test_area():
    # Check that quadrature weights sum to the area of the panel
    mesh = cpt.Sphere().mesh
    for method in [quadpy.c2.schemes['sommariva_05'](), quadpy.c2.schemes["stroud_c2_7_4"]()]:
        mesh.compute_quadrature(method)
        for i_face in range(mesh.nb_faces):
            assert np.isclose(np.sum(mesh.quadrature_points[1][i_face, :]), mesh.faces_areas[i_face], rtol=1e-2)


@pytest.mark.skipif(quadpy is None, reason="quadpy is not installed")
def test_resolution():
    cylinder = cpt.HorizontalCylinder(
        length=5.0, radius=1.0,
        center=(0, 0, -2),
        nr=2, nx=10, ntheta=5,
    )
    # cylinder.show()
    cylinder.add_translation_dof(name="Heave")

    test_matrix = xr.Dataset(coords={
        "omega": np.linspace(0.5, 3.0, 2),
        "radiating_dof": ["Heave"],
    })

    solver = cpt.BEMSolver()

    cylinder.mesh.compute_quadrature(quadpy.c2.schemes['sommariva_01']())
    data_1 = solver.fill_dataset(test_matrix, [cylinder], mesh=True)

    cylinder.mesh.compute_quadrature(quadpy.c2.schemes['sommariva_03']())
    data_3 = solver.fill_dataset(test_matrix, [cylinder], mesh=True)

    assert data_1['quadrature_method'] == "Sommariva 1"
    assert data_3['quadrature_method'] == "Sommariva 3"
    assert np.allclose(data_1["added_mass"].data, data_3["added_mass"].data, rtol=1e-2)
