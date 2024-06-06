import os
import logging

import pytest

import numpy as np
import xarray as xr

import capytaine as cpt

from capytaine.bem.problems_and_results import LinearPotentialFlowProblem, \
        LinearPotentialFlowResult, DiffractionResult, RadiationResult
from capytaine.io.xarray import problems_from_dataset, assemble_dataset
from capytaine.io.legacy import import_cal_file


@pytest.fixture
def sphere():
    sphere = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, -2), radius=1.0, resolution=(10, 20)),
            name="sphere",
            )
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return sphere


@pytest.fixture
def solver():
    solver = cpt.BEMSolver()
    return solver


def test_LinearPotentialFlowProblem(sphere):
    # Without a body
    pb = LinearPotentialFlowProblem(omega=1.0)
    assert pb.omega == 1.0
    assert pb.period == 2*np.pi
    assert pb.wavenumber == 1.0/9.81
    assert pb.wavelength == 9.81*2*np.pi

    assert LinearPotentialFlowProblem(free_surface=np.inf, water_depth=np.inf).water_depth == np.inf
    assert LinearPotentialFlowProblem(free_surface=0.0, water_depth=np.inf).water_depth == np.inf

    pb = LinearPotentialFlowProblem(free_surface=0.0, water_depth=1.0, omega=1.0)
    assert pb.water_depth == 1.0
    assert np.isclose(pb.omega**2, pb.g*pb.wavenumber*np.tanh(pb.wavenumber*pb.water_depth))

    with pytest.raises(NotImplementedError):
        LinearPotentialFlowProblem(free_surface=2.0)

    with pytest.raises(NotImplementedError):
        LinearPotentialFlowProblem(free_surface=np.inf, water_depth=2.0)

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(free_surface=0.0, water_depth=-1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(radiating_dof="Heave")

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(body=cpt.FloatingBody(mesh=cpt.Mesh([], [])))

    # With a body
    pb = LinearPotentialFlowProblem(body=sphere, omega=1.0)
    pb.boundary_condition = sphere.mesh.faces_normals @ (1, 1, 1)
    assert list(pb.influenced_dofs.keys()) == ['Heave']

    pb2 = LinearPotentialFlowProblem(body=sphere, omega=2.0)
    pb2.boundary_condition = sphere.mesh.faces_normals @ (1, 1, 1)
    assert pb < pb2

    # Test transformation to result class
    res = pb.make_results_container()
    assert isinstance(res, LinearPotentialFlowResult)
    assert res.problem is pb
    assert res.omega == pb.omega
    assert res.period == pb.period
    assert res.body is pb.body


def test_mesh_above_the_free_surface(caplog):
    with caplog.at_level(logging.WARNING):
        body = cpt.FloatingBody(mesh=cpt.mesh_sphere(), dofs=cpt.rigid_body_dofs())
        pb = cpt.RadiationProblem(body=body, omega=1.0)
    assert '50 panels above the free surface' in caplog.text
    assert np.all(pb.body.mesh.vertices[:, 2].max() <= 0.0)


def test_mesh_below_the_sea_bottom(caplog):
    with caplog.at_level(logging.WARNING):
        body = cpt.FloatingBody(
                mesh=cpt.mesh_sphere(radius=0.5, center=(0, 0, -1.0)),
                dofs=cpt.rigid_body_dofs()
                )
        pb = cpt.RadiationProblem(body=body, omega=1.0, water_depth=1.0)
    assert '50 panels below the sea bottom' in caplog.text
    assert np.all(pb.body.mesh.vertices[:, 2].min() >= -1.0)


def test_backward_compatibility_with_sea_bottom_argument(caplog):
    with caplog.at_level(logging.WARNING):
        pb = cpt.DiffractionProblem(sea_bottom=-10.0)
        assert pb.water_depth == 10.0
    assert 'water_depth' in caplog.text


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_wavelength(water_depth):
    λ = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(wavelength=λ, water_depth=water_depth).wavelength, λ)


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_wavenumber(water_depth):
    k = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(wavenumber=k, water_depth=water_depth).wavenumber, k)


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_period(water_depth):
    T = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(period=T, water_depth=water_depth).period, T)


def test_setting_too_many_frequencies():
    with pytest.raises(ValueError, match="at most one"):
        cpt.DiffractionProblem(omega=1.0, wavelength=1.0)


def test_diffraction_problem(sphere):
    assert cpt.DiffractionProblem(omega=1.0).body is None

    pb = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, wavenumber=1.0)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    with pytest.raises(TypeError):
        cpt.DiffractionProblem(boundary_conditions=[0, 0, 0])

    assert "DiffractionProblem" in str(cpt.DiffractionProblem(wavenumber=1.0, g=10, rho=1025, free_surface=np.inf))

    res = pb.make_results_container()
    assert isinstance(res, DiffractionResult)


def test_wave_direction_radians_warning(sphere, caplog):
    with caplog.at_level(logging.WARNING):
        cpt.DiffractionProblem(body=sphere, omega=1.0, wave_direction=180)
    assert 'in radians and not in degrees' in caplog.text


def test_radiation_problem(sphere, caplog):
    pb = cpt.RadiationProblem(body=sphere, omega=1.0)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    pb2 = cpt.RadiationProblem(body=sphere, radiating_dof="Heave")
    assert np.all(pb.boundary_condition == pb2.boundary_condition)

    assert "RadiationProblem" in str(cpt.RadiationProblem(g=10, rho=1025, free_surface=np.inf))

    res = pb.make_results_container()
    assert isinstance(res, RadiationResult)
    assert res.added_mass == res.added_masses == {}
    assert res.radiation_damping == res.radiation_dampings == {}

    res = pb.make_results_container(forces={"Heave": 1.0 + 2.0j})
    assert res.added_mass == res.added_masses == {"Heave": 1.0}
    assert res.radiation_damping == res.radiation_dampings == {"Heave": 2.0}


def test_radiation_problem_with_wrong_dof_name():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    with pytest.raises(ValueError, match="'Flurp'"):
        cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Flurp")


@pytest.mark.parametrize("cal_file", ["Nemoh.cal", "Nemoh_v3.cal"])
def test_import_cal_file(cal_file):
    """Test the importation of legacy Nemoh.cal files."""
    current_file_path = os.path.dirname(os.path.abspath(__file__))

    # Non symmetrical body
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "NonSymmetrical", cal_file)
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*41+41
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.water_depth == np.inf
        assert isinstance(problem.body, cpt.FloatingBody)
        assert problem.body.nb_dofs == 6
        assert problem.body.mesh.nb_vertices == 299  # Duplicate vertices are removed during import.
        assert problem.body.mesh.nb_faces == 280
        assert problem.omega in np.linspace(0.1, 2.0, 41)
        if isinstance(problem, cpt.DiffractionProblem):
            assert problem.wave_direction == 0.0

    # Symmetrical cylinder
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "Cylinder", cal_file)
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*2+2
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.water_depth == np.inf
        assert isinstance(problem.body.mesh, cpt.ReflectionSymmetricMesh)
        assert isinstance(problem.body.mesh[0], cpt.Mesh)
        assert problem.body.nb_dofs == 6
        # assert problem.body.mesh.nb_vertices == 2*540
        assert problem.body.mesh.nb_faces == 2*300
        assert problem.omega == 0.1 or problem.omega == 2.0
        if isinstance(problem, cpt.DiffractionProblem):
            assert problem.wave_direction == 0.0


def test_results():
    assert isinstance(LinearPotentialFlowProblem().make_results_container(), LinearPotentialFlowResult)

    pb = cpt.DiffractionProblem(g=10, rho=1023, free_surface=np.inf)
    res = DiffractionResult(pb)
    assert res.g == pb.g == 10
    assert "DiffractionResult" in str(res)

    pb = cpt.RadiationProblem(g=10, rho=1023, free_surface=np.inf)
    res = RadiationResult(pb)
    assert res.g == pb.g == 10
    assert "RadiationResult" in str(res)


def test_problems_from_dataset(sphere):
    dset = xr.Dataset(coords={'omega': [0.5, 1.0, 1.5],
                              'radiating_dof': ["Heave"],
                              'body_name': ["sphere"],
                              'wave_direction': [0.0],
                              'water_depth': [np.inf]})

    problems = problems_from_dataset(dset, [sphere])
    assert cpt.RadiationProblem(body=sphere, omega=0.5, radiating_dof="Heave") in problems
    assert len(problems) == 6
    assert len([problem for problem in problems if isinstance(problem, cpt.DiffractionProblem)]) == 3

    dset = xr.Dataset(coords={'omega': [0.5, 1.0, 1.5],
                              'wave_direction': [0.0],
                              'body_name': ["cube"]})
    with pytest.raises(AssertionError):
        problems_from_dataset(dset, [sphere])

    shifted_sphere = sphere.translated_y(5.0, name="shifted_sphere")
    dset = xr.Dataset(coords={'omega': [0.5, 1.0, 1.5],
                              'radiating_dof': ["Heave"],
                              'wave_direction': [0.0]})
    problems = problems_from_dataset(dset, [sphere, shifted_sphere])
    assert cpt.RadiationProblem(body=sphere, omega=0.5, radiating_dof="Heave") in problems
    assert cpt.RadiationProblem(body=shifted_sphere, omega=0.5, radiating_dof="Heave") in problems
    assert len(problems) == 12


def test_problems_from_dataset_with_wavelength(sphere):
    dset = xr.Dataset(coords={'wavelength': [12.0], 'radiating_dof': ["Heave"]})
    problems = problems_from_dataset(dset, sphere)
    for pb in problems:
        assert np.isclose(pb.wavelength, 12.0)


def test_problems_from_dataset_with_too_many_frequencies(sphere):
    dset = xr.Dataset(coords={'wavelength': [12.0], 'period': [3.0], 'radiating_dof': ["Heave"]})
    with pytest.raises(ValueError, match="at most one"):
        problems_from_dataset(dset, sphere)


def test_problems_from_dataset_without_list(sphere):
    dset = xr.Dataset(coords={'omega': 1.5, 'radiating_dof': "Heave"})
    problems = problems_from_dataset(dset, sphere)
    assert all(pb.omega == 1.5 for pb in problems)
    assert all(pb.radiating_dof == "Heave" for pb in problems)


def test_problems_from_dataset_without_list_with_too_many_frequencies(sphere):
    dset = xr.Dataset(coords={'omega': 1.5, 'period': 1.5, 'radiating_dof': "Heave"})
    with pytest.raises(ValueError, match="at most one"):
        problems_from_dataset(dset, sphere)


@pytest.mark.parametrize("method", ['indirect','direct'])
def test_assemble_dataset(sphere, solver, method):
    pb_1 = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0)
    res_1 = solver.solve(pb_1, method=method)
    ds1 = assemble_dataset([res_1])
    assert "Froude_Krylov_force" in ds1

    pb_2 = cpt.RadiationProblem(body=sphere, radiating_dof="Heave", omega=1.0)
    res_2 = solver.solve(pb_2, method=method)
    ds2 = assemble_dataset([res_2])
    assert "added_mass" in ds2

    ds12 = assemble_dataset([res_1, res_2])
    assert "Froude_Krylov_force" in ds12
    assert "added_mass" in ds12


def test_fill_dataset(sphere, solver):
    sphere.add_all_rigid_body_dofs()
    test_matrix = xr.Dataset(coords={'omega': [1.0, 2.0, 3.0], 'wave_direction': [0, np.pi/2], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, [sphere])
    assert dataset['added_mass'].data.shape == (3, 1, 6)
    assert dataset['Froude_Krylov_force'].data.shape == (3, 2, 6)


def test_fill_dataset_with_wavenumbers(sphere, solver):
    k_range = np.linspace(1.0, 3.0, 3)
    test_matrix = xr.Dataset(coords={'wavenumber': k_range, 'wave_direction': [0, np.pi/2], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, [sphere])
    np.testing.assert_allclose(dataset.coords['wavenumber'], k_range)
    assert set(dataset.added_mass.dims) == {'wavenumber', 'radiating_dof', 'influenced_dof'}
    assert set(dataset.wavenumber.dims) == {'wavenumber'}
    assert set(dataset.wavelength.dims) == {'wavenumber'}
    assert set(dataset.omega.dims)      == {'wavenumber'}
    assert set(dataset.period.dims)     == {'wavenumber'}


def test_fill_dataset_with_periods(sphere, solver):
    T_range = np.linspace(1.0, 3.0, 3)
    test_matrix = xr.Dataset(coords={'period': T_range, 'wave_direction': [0, np.pi/2], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, [sphere])

    np.testing.assert_allclose(sorted(dataset.coords['period']), sorted(T_range))
    assert set(dataset.added_mass.dims) == {'period', 'radiating_dof', 'influenced_dof'}
    assert set(dataset.wavenumber.dims) == {'period'}
    assert set(dataset.wavelength.dims) == {'period'}
    assert set(dataset.omega.dims)      == {'period'}
    assert set(dataset.period.dims)     == {'period'}


def test_fill_dataset_with_wavenumbers_and_several_water_depths(sphere, solver):
    k_range = np.linspace(1.0, 3.0, 3)
    test_matrix = xr.Dataset(coords={
        'wavenumber': k_range, 'radiating_dof': ['Heave'], 'water_depth': [4.0, 6.0],
    })
    dataset = solver.fill_dataset(test_matrix, [sphere])

    np.testing.assert_allclose(dataset.coords['wavenumber'], k_range)
    assert set(dataset.added_mass.dims) == {'wavenumber', 'radiating_dof', 'influenced_dof', 'water_depth'}
    assert set(dataset.wavenumber.dims) == {'wavenumber'}
    assert set(dataset.wavelength.dims) == {'wavenumber'}
    assert set(dataset.omega.dims)      == {'wavenumber', 'water_depth'}
    assert set(dataset.period.dims)     == {'wavenumber', 'water_depth'}
