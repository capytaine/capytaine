import pytest

import numpy as np
import xarray as xr

import capytaine as cpt


@pytest.fixture
def sphere():
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body


def test_limit_frequencies(sphere):
    # TODO: do for direct and indirect
    """Test if how the solver answers when asked for frequency of 0 or ∞."""
    solver = cpt.BEMSolver()

    solver.solve(cpt.RadiationProblem(body=sphere, omega=0.0, water_depth=np.inf))

    with pytest.raises(NotImplementedError):
        solver.solve(cpt.RadiationProblem(body=sphere, omega=0.0, water_depth=1.0))

    solver.solve(cpt.RadiationProblem(body=sphere, omega=np.inf, water_depth=np.inf))

    with pytest.raises(NotImplementedError):
        solver.solve(cpt.RadiationProblem(body=sphere, omega=np.inf, water_depth=10))


def test_radiation_damping_value(sphere):
    solver = cpt.BEMSolver()
    res = solver.solve(cpt.RadiationProblem(body=sphere, omega=0.0, water_depth=np.inf))
    assert res.radiation_damping['Surge'] == 0.0
    res = solver.solve(cpt.RadiationProblem(body=sphere, omega=np.inf, water_depth=np.inf))
    assert res.radiation_damping['Surge'] == 0.0


def test_limit_frequencies_with_symmetries():
    mesh = cpt.mesh_parallelepiped(reflection_symmetry=True).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Surge")
    pb = cpt.RadiationProblem(body=body, omega=0.0)
    solver = cpt.BEMSolver()
    res = solver.solve(pb, keep_details=True)
    assert isinstance(res.added_mass['Surge'], float)


def test_zero_frequency_datasets(sphere):
    test_matrix = xr.Dataset(coords={
        "wavenumber": np.linspace(0.0, 0.1, 3),
        "radiating_dof": list(sphere.dofs),
        })
    solver = cpt.BEMSolver()
    solver.fill_dataset(test_matrix, sphere)


def test_infinite_frequency_datasets(sphere):
    test_matrix = xr.Dataset(coords={
        "wavelength": np.linspace(0.0, 0.1, 3),
        "radiating_dof": list(sphere.dofs),
        })
    solver = cpt.BEMSolver()
    solver.fill_dataset(test_matrix, sphere)


def test_no_warning_mesh_resolution_at_zero_wavelength(sphere, caplog):
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=sphere, wavelength=0)
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "resolution " not in caplog.text


def test_dataset_with_zero_frequency_including_radiation_and_diffraction(sphere):
    test_matrix = xr.Dataset(coords={
        'omega': [0.0, 1.0, np.inf],
        'wave_direction': [0.0],
        'radiating_dof': list(sphere.dofs),
    })
    solver = cpt.BEMSolver()
    ds = solver.fill_dataset(test_matrix, sphere)
    assert np.all(np.isnan(ds.diffraction_force.sel(omega=0.0)))
    assert not np.any(np.isnan(ds.diffraction_force.sel(omega=1.0)))
    assert np.all(np.isnan(ds.diffraction_force.sel(omega=np.inf)))
    assert not np.any(np.isnan(ds.added_mass))
    assert not np.any(np.isnan(ds.radiation_damping))


def test_dataset_with_zero_frequency_diffraction_only(sphere):
    test_matrix = xr.Dataset(coords={
        'omega': [0.0, 1.0, np.inf],
        'wave_direction': [0.0, np.pi],
    })
    solver = cpt.BEMSolver()
    with pytest.raises(ValueError):
        solver.fill_dataset(test_matrix, sphere)
