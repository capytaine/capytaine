import pytest
import numpy as np
import xarray as xr

import capytaine as cpt
from capytaine.post_pro import irf

@pytest.fixture
def sphere_fb():
    m = 1.866e+03
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(3, 12)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), mass=m)
    body.inertia_matrix = body.add_dofs_labels_to_matrix(
            [[ m ,  0.0, 0.0, 0.0,        0.0,       0.0       ],
             [ 0.0, m,   0.0, 0.0,        0.0,       0.0       ],
             [ 0.0, 0.0, m,   0.0,        0.0,       0.0       ],
             [ 0.0, 0.0, 0.0, 4.469e+02,  9.676e-31, -2.757e-14],
             [ 0.0, 0.0, 0.0, 9.676e-31,  4.469e+02, 3.645e-15 ],
             [ 0.0, 0.0, 0.0, -2.757e-14, 3.645e-15, 6.816e+02 ]]
            )
    body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(
            [[    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
             [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
             [    0.   ,     0.   , 29430.   ,     0.   ,     0.   ,     0.   ],
             [    0.   ,     0.   ,     0.   ,   328.573,     0.   ,     0.   ],
             [    0.   ,     0.   ,     0.   ,     0.   ,   328.573,     0.   ],
             [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ]]
            )
    return body


def test_rao_sphere_all(sphere_fb):
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.5, 10.0, 5),
        'wave_direction': [0],
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, sphere_fb,
                               hydrostatics=True, mesh=True,
                               wavelength=True, wavenumber=True)

    RAO = cpt.post_pro.rao(data)

    time_series = np.arange(0,50,0.1)
    kr = irf.rad_kernel(data, time_series)

    assert kr.shape == (len(time_series), 
                        data.radiating_dof.size,
                        data.influenced_dof.size)
    assert RAO.radiating_dof.size == 6
    assert data.inertia_matrix.shape == (6,6)
    assert np.all(data.inertia_matrix.values == sphere_fb.inertia_matrix.values)
    assert data.hydrostatic_stiffness.shape == (6,6)
    assert np.all(data.hydrostatic_stiffness.values == sphere_fb.hydrostatic_stiffness.values)
    # # assert RAO == ? # TODO could test against known results


def test_rao_from_wavelengths(sphere_fb):
    # From https://github.com/capytaine/capytaine/issues/316
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
        'wavelength': np.linspace(0.5, 10.0, 5),
        'wave_direction': [0],
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, sphere_fb, hydrostatics=True)
    Zi = cpt.post_pro.impedance(data)
    RAO = cpt.post_pro.rao(data)


@pytest.fixture
def sphere_heave_data(sphere_fb):
    sphere_fb.keep_only_dofs(['Heave'])

    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
          'omega': np.linspace(0.5, 10.0, 5),
          'wave_direction': [0],
          'radiating_dof': list(sphere_fb.dofs.keys()),
          })

    data = solver.fill_dataset(test_matrix, [sphere_fb], hydrostatics=True)

    return data


def test_impedance_sphere_heave(sphere_heave_data):
    Zi = cpt.post_pro.impedance(sphere_heave_data)


def test_malformed_dataset(sphere_heave_data):
    data = sphere_heave_data.drop_vars("wavelength")
    data = data.expand_dims({"wavelength": 1})
    # data has both an "omega" dimension and a "wavelength" dimension
    with pytest.raises(ValueError):
        RAO = cpt.post_pro.rao(data)


def test_rao_sphere_heave_indirect(sphere_heave_data):
    RAO = cpt.post_pro.rao(sphere_heave_data)
    assert RAO.radiating_dof.size == 1
    assert sphere_heave_data.inertia_matrix.size == 1
    assert sphere_heave_data.hydrostatic_stiffness.size == 1
    # assert RAO == ? # TODO could test against known results
