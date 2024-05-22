import pytest
import numpy as np
from scipy.special import j0
import capytaine as cpt


def test_infinite_depth_gf():
    gf = cpt.Delhommeau()
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                args[3])
    xi = np.array([0.0, 0.0, -0.5])
    xj = np.array([0.0, 1.0, -0.5])
    k = 2.0
    g_lf, nablag_lf = wave_part(xi, xj, k, 0)
    g_hf, nablag_hf = wave_part(xi, xj, k, 1)
    r = k * np.linalg.norm(xi[:2] - xj[:2])
    z = k * (xi[2] + xj[2])
    assert g_hf.real + 2*k/np.hypot(r, z) == pytest.approx(g_lf.real, rel=1e-4)
    assert g_lf.imag == pytest.approx(g_hf.imag)
    assert g_lf.imag == pytest.approx(2*np.pi*k*np.exp(z)*j0(r), rel=1e-4)


@pytest.mark.parametrize("gf_singularities", ["low_freq", "high_freq"])
def test_infinite_depth_gf_finite_differences(gf_singularities):
    gf = cpt.Delhommeau(gf_singularities=gf_singularities)
    def wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                gf.gf_singularities_index)[0]
    def nabla_wave_part(*args):
        return gf.fortran_core.green_wave.wave_part_infinite_depth(
                *args[:3], gf.tabulation_nb_integration_points, gf.tabulation_grid_shape_index,
                gf.tabulated_r_range, gf.tabulated_z_range, gf.tabulated_integrals,
                gf.gf_singularities_index)[1]
    k = 1.0
    xi = np.array([0.0, 0.0, -0.5])
    xj = np.array([1.0, 1.0, -0.5])
    x_eps = np.array([1e-3, 0.0, 0.0])
    y_eps = np.array([0.0, 1e-3, 0.0])
    z_eps = np.array([0.0, 0.0, 1e-3])
    g = wave_part(xi, xj, k)
    g_x = (wave_part(xi + x_eps, xj, k) - g)/np.linalg.norm(x_eps)
    g_y = (wave_part(xi + y_eps, xj, k) - g)/np.linalg.norm(y_eps)
    g_z = (wave_part(xi + z_eps, xj, k) - g)/np.linalg.norm(z_eps)
    nabla_g = np.array([g_x, g_y, g_z])
    nablag_ref = nabla_wave_part(xi, xj, k)
    assert nablag_ref == pytest.approx(nabla_g, abs=1e-2)
