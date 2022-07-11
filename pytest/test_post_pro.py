import pytest
import numpy as np
import xarray as xr
from capytaine import BEMSolver
from capytaine.bodies.predefined.spheres import Sphere
from capytaine.post_pro import rao
from capytaine.post_pro import impedance

f = np.linspace(0.1, 2.0)
omega = 2*np.pi*f
rho_water = 1e3
r = 1
m = 1.866e+03

M = np.array([
       [ m,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00],
       [ 0.000e+00,  m,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00],
       [ 0.000e+00,  0.000e+00,  m,  0.000e+00,  0.000e+00,  0.000e+00],
       [ 0.000e+00,  0.000e+00,  0.000e+00,  4.469e+02,  9.676e-31, -2.757e-14],
       [ 0.000e+00,  0.000e+00,  0.000e+00,  9.676e-31,  4.469e+02,  3.645e-15],
       [ 0.000e+00,  0.000e+00,  0.000e+00, -2.757e-14,  3.645e-15,  6.816e+02]])

kHS = np.array([
    [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
    [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
    [    0.   ,     0.   , 29430.   ,     0.   ,     0.   ,     0.   ],
    [    0.   ,     0.   ,     0.   ,   328.573,     0.   ,     0.   ],
    [    0.   ,     0.   ,     0.   ,     0.   ,   328.573,     0.   ],
    [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ]])


@pytest.fixture
def sphere_fb():
    sphere = Sphere(radius=r, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_all_rigid_body_dofs()

    sphere.inertia_matrix = sphere.add_dofs_labels_to_matrix(M)
    sphere.hydrostatic_stiffness = sphere.add_dofs_labels_to_matrix(kHS)

    return sphere

def test_rao_sphere_all(sphere_fb):

    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
        'rho': rho_water,
        'water_depth': [np.infty],
        'omega': omega,
        'wave_direction': 0,
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, [sphere_fb],
                               hydrostatics=True,
                               mesh=True,
                               wavelength=True,
                               wavenumber=True)

    RAO = rao(data)

    assert RAO.radiating_dof.size == 6
    assert data.inertia_matrix.shape == (6,6)
    assert np.all(data.inertia_matrix.values == M)
    assert data.hydrostatic_stiffness.shape == (6,6)
    assert np.all(data.hydrostatic_stiffness.values == kHS)
    # # assert RAO == ? # TODO could test against known results

@pytest.fixture
def sphere_heave_data(sphere_fb):
    sphere_fb.keep_only_dofs(['Heave'])

    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
          'rho': rho_water,
          'water_depth': [np.infty],
          'omega': omega,
          'wave_direction': 0,
          'radiating_dof': list(sphere_fb.dofs.keys()),
          })

    data = solver.fill_dataset(test_matrix, [sphere_fb],
                                 hydrostatics=True,
                                 mesh=True,
                                 wavelength=True,
                                 wavenumber=True)

    return data

def test_impedance_sphere_heave(sphere_heave_data):
    Zi = impedance(sphere_heave_data)


def test_rao_sphere_heave_indirect(sphere_heave_data):
    RAO = rao(sphere_heave_data)
    assert RAO.radiating_dof.size == 1
    assert sphere_heave_data.inertia_matrix.size == 1
    assert sphere_heave_data.inertia_matrix.values == m
    assert sphere_heave_data.hydrostatic_stiffness.size == 1
    assert sphere_heave_data.hydrostatic_stiffness.values == kHS[2,2]
    # assert RAO == ? # TODO could test against known results
