import pytest
import numpy as np
import xarray as xr
import capytaine as cpt


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

@pytest.fixture
def solver():
    return cpt.BEMSolver()

def test_rao_sphere_all(sphere_fb, solver):
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.5, 10.0, 5),
        'wave_direction': [0],
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, sphere_fb,
                               hydrostatics=True, mesh=True)

    RAO = cpt.post_pro.rao(data)

    assert RAO.radiating_dof.size == 6
    assert data.inertia_matrix.shape == (6,6)
    assert np.all(data.inertia_matrix.values == sphere_fb.inertia_matrix.values)
    assert data.hydrostatic_stiffness.shape == (6,6)
    assert np.all(data.hydrostatic_stiffness.values == sphere_fb.hydrostatic_stiffness.values)

    Fex_ = xr.dot(cpt.post_pro.rao_transfer_function(data), RAO, dims="radiating_dof")
    assert np.allclose(Fex_.values.squeeze(), data["excitation_force"].values.squeeze())
    # # assert RAO == ? # TODO could test against known results


def test_rao_from_wavelengths(sphere_fb, solver):
    # From https://github.com/capytaine/capytaine/issues/316
    test_matrix = xr.Dataset(coords={
        'wavelength': np.linspace(0.5, 10.0, 5),
        'wave_direction': [0],
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, sphere_fb, hydrostatics=True)
    Zi = cpt.post_pro.impedance(data)
    RAO = cpt.post_pro.rao(data)


def test_rao_several_water_depth(sphere_fb, solver):
    # From https://github.com/capytaine/capytaine/issues/405
    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.5, 10.0, 3),
        'wave_direction': [0],
        'water_depth': [np.inf, 10.0],
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })

    data = solver.fill_dataset(test_matrix, sphere_fb, hydrostatics=True)
    rao = cpt.post_pro.rao(data)
    assert "water_depth" in rao.dims


@pytest.fixture
def sphere_heave_data(solver, sphere_fb):
    body = sphere_fb.with_only_dofs(['Heave'])
    test_matrix = xr.Dataset(coords={
          'omega': np.linspace(0.5, 10.0, 5),
          'wave_direction': [0],
          'radiating_dof': list(body.dofs.keys()),
          })

    data = solver.fill_dataset(test_matrix, [body], hydrostatics=True)

    return data


def test_impedance_sphere_heave(sphere_heave_data):
    Zi = cpt.post_pro.impedance(sphere_heave_data)


def test_rao_sphere_heave_indirect(sphere_heave_data):
    RAO = cpt.post_pro.rao(sphere_heave_data)
    assert RAO.radiating_dof.size == 1
    assert sphere_heave_data.inertia_matrix.size == 1
    assert sphere_heave_data.hydrostatic_stiffness.size == 1
    # assert RAO == ? # TODO could test against known results


def test_asymmetric_matrices():
    omega = 2.0
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=["Surge", "Heave"]))
    M = mesh.disp_mass() * np.array(
            [[ 1.0 , 0.0 ],
             [ 0.0 , 1.0 ]]
            )
    body.inertia_matrix = body.add_dofs_labels_to_matrix(M)
    S = mesh.disp_mass() * np.array(
            [[ 1.0 , 10.0 ],
             [ 0.0 , 1.0 ]]
            )
    body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(S)
    test_matrix = xr.Dataset(coords={
        'omega': [omega],
        'wave_direction': [0],
        'radiating_dof': list(body.dofs),
        })
    solver = cpt.BEMSolver()
    data = solver.fill_dataset(test_matrix, body, hydrostatics=True)
    X = cpt.post_pro.rao(data)
    assert np.allclose(
        (cpt.post_pro.rao_transfer_function(data) @ X).values.squeeze(),
        data["excitation_force"].values.squeeze()
    )
    A = data["added_mass"].sel(omega=omega).values
    B = data["radiation_damping"].sel(omega=omega).values
    Fex = data["excitation_force"].sel(omega=omega, wave_direction=0.0).values
    H = -omega**2*(M + A) - 1j * omega * B + S
    X_ref = np.linalg.solve(H, Fex)
    assert np.allclose(X_ref, X.values.squeeze())
