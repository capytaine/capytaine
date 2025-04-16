import pytest

import numpy as np
import xarray as xr

import capytaine as cpt

from capytaine.io.xarray import problems_from_dataset


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


#######################################################################
#                       Problems from datasets                        #
#######################################################################

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


#######################################################################
#                          Assemble matrices                          #
#######################################################################

def test_assemble_matrices(sphere, solver):
    pbs = [cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0),
           cpt.RadiationProblem(body=sphere, wave_direction=1.0, radiating_dof="Heave")]
    res = solver.solve_all(pbs)
    A, B, F = cpt.assemble_matrices(res)
    assert A.shape == (1, 1)
    assert A.dtype == np.float64
    assert B.shape == (1, 1)
    assert B.dtype == np.float64
    assert F.shape == (1,)
    assert F.dtype == np.complex128


def test_assemble_matrices_rad_only(sphere, solver):
    pbs = [cpt.RadiationProblem(body=sphere, wave_direction=1.0, radiating_dof="Heave")]
    res = solver.solve_all(pbs)
    A, B, F = cpt.assemble_matrices(res)
    assert A.shape == (1, 1)
    assert A.dtype == np.float64
    assert B.shape == (1, 1)
    assert B.dtype == np.float64
    assert F is None


def test_assemble_matrices_dif_only(sphere, solver):
    pbs = [cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0)]
    res = solver.solve_all(pbs)
    A, B, F = cpt.assemble_matrices(res)
    assert A is None
    assert B is None
    assert F.shape == (1,)
    assert F.dtype == np.complex128


def test_assemble_matrices_no_data():
    with pytest.raises(ValueError):
        cpt.assemble_matrices([])

#######################################################################
#                          Assemble dataset                           #
#######################################################################

def test_assemble_dataset(sphere, solver):
    pb_1 = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0)
    res_1 = solver.solve(pb_1)
    ds1 = cpt.assemble_dataset([res_1])
    assert "diffraction_force" in ds1
    assert "added_mass" not in ds1

    pb_2 = cpt.RadiationProblem(body=sphere, radiating_dof="Heave", omega=1.0)
    res_2 = solver.solve(pb_2)
    ds2 = cpt.assemble_dataset([res_2])
    assert "added_mass" in ds2
    assert "diffraction_force" not in ds2

    ds12 = cpt.assemble_dataset([res_1, res_2])
    assert "diffraction_force" in ds12
    assert "added_mass" in ds12


def test_assemble_dataset_with_nans(sphere):
    pb = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0)
    res = pb.make_results_container(forces={dof: np.nan for dof in pb.influenced_dofs})
    ds = cpt.assemble_dataset([res])
    assert 'diffraction_force' in ds
    assert np.all(np.isnan(ds.diffraction_force.values))


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


@pytest.fixture
def broken_bem_solver():
    ref_gf = cpt.Delhommeau()
    class BrokenGreenFunction:
        def evaluate(self, m1, m2, fs, wd, wavenumber, *args, **kwargs):
            if wavenumber < 2.0:
                raise NotImplementedError("I'm potato")
            else:
                return ref_gf.evaluate(m1, m2, fs, wd, wavenumber, *args, **kwargs)
    broken_bem_solver = cpt.BEMSolver(green_function=BrokenGreenFunction())
    return broken_bem_solver


def test_failed_resolution_in_dataset(broken_bem_solver, sphere):
    test_matrix = xr.Dataset(coords={"wavenumber": np.linspace(0.1, 5.0, 5), "wave_direction": [0.0], "radiating_dof": ["Heave"]})
    ds = broken_bem_solver.fill_dataset(test_matrix, sphere)
    assert len(ds.wavenumber) == 5
    assert np.any(np.isnan(ds.added_mass.values))
