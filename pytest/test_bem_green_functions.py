"""Tests for the computation of the Green function using Fortran routines."""

import pytest
from pytest import approx

import numpy as np
import capytaine as cpt


gfs = [
        cpt.Delhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        cpt.XieDelhommeau(tabulation_nr=328, tabulation_nz=46, tabulation_nb_integration_points=251, tabulation_grid_shape="legacy"),
        ]


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_infinite_depth(gf):
    k = 1.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xi, xj, k,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals, gf.gf_singularities_index
    )
    g2, dg2 = gf.fortran_core.green_wave.wave_part_infinite_depth(
        xj, xi, k,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals, gf.gf_singularities_index
    )
    assert g1 == approx(g2)
    assert dg1[0:2] == approx(-dg2[0:2])
    assert dg1[2] == approx(dg2[2])


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_finite_depth_no_prony(gf):
    k = 1.0
    depth = 5.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xi, xj, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        np.zeros(1), np.zeros(1), 1
    )
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xj, xi, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        np.zeros(1), np.zeros(1), 1
    )
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


@pytest.mark.parametrize("gf", gfs)
def test_symmetry_of_the_green_function_finite_depth(gf):
    k = 1.0
    depth = 10.0
    xi = np.array([0.0, 0.0, -1.0])
    xj = np.array([1.0, 1.0, -2.0])
    ambda, a, nexp = gf.fortran_core.old_prony_decomposition.lisc(k*depth*np.tanh(k*depth), k*depth)
    g1, dg1_sym, dg1_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xi, xj, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        ambda, a, 31
    )
    g2, dg2_sym, dg2_antisym = gf.fortran_core.green_wave.wave_part_finite_depth(
        xj, xi, k, depth,
        gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
        gf.tabulated_r_range, gf.tabulated_z_range,
        gf.tabulated_integrals,
        ambda, a, 31
    )
    assert g1 == approx(g2)
    assert dg1_sym == approx(dg2_sym)
    assert dg1_antisym == approx(-dg2_antisym)


def test_floating_point_precision():
    assert cpt.Delhommeau(floating_point_precision="float64").tabulated_integrals.dtype == np.float64
    assert cpt.Delhommeau(floating_point_precision="float32").tabulated_integrals.dtype == np.float32


def test_no_tabulation():
    mesh = cpt.Sphere().mesh.keep_immersed_part()
    tabed_gf = cpt.Delhommeau()
    untabed_gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0)
    assert np.allclose(untabed_gf.evaluate(mesh, mesh)[0], tabed_gf.evaluate(mesh, mesh)[0], atol=1e-2)


# def test_panels_near_free_surface():
#     area = 4.0
#     def gf_at_depth(d):
#         mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), resolution=(1, 1), center=(0, 0, d))
#         x = np.array([[0, 0, 0]])
#         S, _ = cpt.Delhommeau().evaluate(x, mesh, free_surface=0, water_depth=np.infty, wavenumber=3.0)
#         return S
#     gf_at_depth(0.0)
#     gf_at_depth(1e-9)
#     gf_at_depth(1e-5)
#     gf_at_depth(1e-2)


@pytest.mark.parametrize("gf", gfs)
def test_exact_integration_of_rankine_terms(gf):
    center = np.array([0.0, 0.0, 0.0])
    area = 1.0

    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(1, 1))
    rankine_g_once = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    mesh = cpt.mesh_rectangle(size=(np.sqrt(area), np.sqrt(area)), center=center, resolution=(11, 11))
    # Use odd number of panels to avoid having the center on a corner of a panel, which is not defined for the strong singularity of the derivative.
    rankine_g_parts = gf.evaluate(center.reshape(1, -1), mesh, wavenumber=0.0)[0]
    assert rankine_g_once == pytest.approx(rankine_g_parts.sum(), abs=1e-2)


def test_exact_integration_with_nb_integration_points():
    # Test that the tabulation_nb_integration_points is actually taken into account.
    mesh = cpt.mesh_rectangle(center=(0, 0, -1), resolution=(1, 1))
    point = np.array([[100, 0, -1]])
    gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nb_integration_points=11)
    val1 = gf.evaluate(point, mesh)
    gf = cpt.Delhommeau(tabulation_nr=0, tabulation_nb_integration_points=101)
    val2 = gf.evaluate(point, mesh)
    assert np.abs(val1[0] - val2[0]) > 1e-1
