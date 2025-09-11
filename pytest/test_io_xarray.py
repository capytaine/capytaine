import pytest

import numpy as np
import xarray as xr

import capytaine as cpt

from capytaine.io.xarray import problems_from_dataset, separate_complex_values, merge_complex_values


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


def test_problems_from_dataset_incomplete_test_matrix():
    body = cpt.FloatingBody(cpt.mesh_sphere())
    test_matrix = xr.Dataset(coords={
        "omega": np.linspace(0, 1, 2),
        })
    with pytest.raises(ValueError):
        problems_from_dataset(test_matrix, body)


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
#                         Assemble dataframes                         #
#######################################################################

def test_assemble_dataframe(sphere, solver):
    pb_1 = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, omega=1.0)
    res_1 = solver.solve(pb_1)
    df1 = cpt.assemble_dataframe([res_1])
    assert "diffraction_force" in df1
    assert "added_mass" not in df1

    pb_2 = cpt.RadiationProblem(body=sphere, radiating_dof="Heave", omega=1.0)
    res_2 = solver.solve(pb_2)
    df2 = cpt.assemble_dataframe([res_2])
    assert "added_mass" in df2
    assert "diffraction_force" not in df2

    df12 = cpt.assemble_dataframe([res_1, res_2])
    assert "diffraction_force" in df12
    assert "added_mass" in df12


def test_assemble_dataframe_with_infinite_free_surface(sphere, solver):
    pb = cpt.RadiationProblem(body=sphere, free_surface=np.inf)
    res = solver.solve(pb)
    df = cpt.assemble_dataframe([res])
    assert "free_surface" in df


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


def test_xarray_dataset_with_more_data():
    # Store some mesh data when several bodies in dataset
    bodies = [
        cpt.FloatingBody(cpt.mesh_sphere(radius=1, resolution=(3, 3)).immersed_part(), name="sphere_1"),
        cpt.FloatingBody(cpt.mesh_sphere(radius=2, resolution=(5, 3)).immersed_part(), name="sphere_2"),
        cpt.FloatingBody(cpt.mesh_sphere(radius=3, resolution=(7, 3)).immersed_part(), name="sphere_3"),
    ]
    for body in bodies:
        body.add_translation_dof(name="Heave")

    problems = [cpt.RadiationProblem(body=b, radiating_dof="Heave", omega=1.0) for b in bodies]
    results = cpt.BEMSolver().solve_all(problems)

    ds = cpt.assemble_dataset(results, mesh=True)
    assert 'nb_faces' in ds.coords
    assert set(ds.coords['nb_faces'].values) == set([b.mesh.nb_faces for b in bodies])


def test_variables_attrs(sphere, solver):
    pb = cpt.RadiationProblem(body=sphere, omega=1.0, radiating_dof="Heave")
    ds = cpt.assemble_dataset([solver.solve(pb)])
    assert 'long_name' in ds.omega.attrs
    assert 'long_name' in ds.wavenumber.attrs
    assert 'long_name' in ds.freq.attrs
    assert 'long_name' in ds.added_mass.attrs


def test_assemble_dataset_with_infinite_free_surface(caplog, sphere, solver):
    pb = cpt.RadiationProblem(body=sphere, free_surface=np.inf)
    res = solver.solve(pb)
    with caplog.at_level("WARNING"):
        ds = cpt.assemble_dataset([res])
    assert "ignored" in caplog.text
    assert len(ds) == 0


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


def test_fill_dataset_with_freqs(sphere, solver):
    f_range = np.linspace(0.1, 1, 3)
    test_matrix = xr.Dataset(coords={'freq': f_range, 'wave_direction': [0, np.pi/2], 'radiating_dof': ['Heave']})
    dataset = solver.fill_dataset(test_matrix, [sphere])
    np.testing.assert_allclose(dataset.coords['freq'], f_range)
    assert set(dataset.added_mass.dims) == {'freq', 'radiating_dof', 'influenced_dof'}
    assert set(dataset.freq.dims) == {'freq'}
    assert set(dataset.wavelength.dims) == {'freq'}
    assert set(dataset.omega.dims)      == {'freq'}
    assert set(dataset.period.dims)     == {'freq'}


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
        exportable_settings = {}
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


#######################################################################
#                       separate_complex_values                       #
#######################################################################

def test_remove_complex_values():
    original_dataset = xr.Dataset(
        data_vars={
            'variable1': (['x', 'y'], np.random.rand(4, 5)),
            'variable2': (['x', 'y', 'z'], np.random.rand(4, 5, 3) + 1j * np.random.rand(4, 5, 3)),
            'variable3': (['y', 'z'], np.random.rand(5, 3) + 1j * np.random.rand(5, 3)),
        },
        coords={
            'x': (['x'], np.linspace(0, 10, 4))
        })

    real_dataset = separate_complex_values(original_dataset)
    assert np.allclose(real_dataset['variable3'].sel(complex='re').data,
                       np.real(original_dataset['variable3'].data),
                       atol=0.1)
    assert set(original_dataset.dims) == set(real_dataset.dims) - {'complex'}

    complex_dataset = merge_complex_values(real_dataset)
    assert set(original_dataset.dims) == set(complex_dataset.dims)
