
import xarray as xr
import numpy as np
import capytaine as cpt
import pytest
import os

from capytaine.io.xarray import problems_from_dataset, separate_complex_values, merge_complex_values


def test_incomplete_test_matrix():
    body = cpt.Sphere()
    test_matrix = xr.Dataset(coords={
        "omega": np.linspace(0, 1, 2),
        })
    with pytest.raises(ValueError):
        problems_from_dataset(test_matrix, body)


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


def test_xarray_dataset_with_more_data():
    # Store some mesh data when several bodies in dataset
    bodies = [
        cpt.Sphere(radius=1, ntheta=3, nphi=3, name="sphere_1"),
        cpt.Sphere(radius=2, ntheta=5, nphi=3, name="sphere_2"),
        cpt.Sphere(radius=3, ntheta=7, nphi=3, name="sphere_3"),
    ]
    for body in bodies:
        body.keep_immersed_part()
        body.add_translation_dof(name="Heave")

    problems = [cpt.RadiationProblem(body=b, radiating_dof="Heave", omega=1.0) for b in bodies]
    results = cpt.BEMSolver().solve_all(problems)

    ds = cpt.assemble_dataset(results, mesh=True)
    assert 'nb_faces' in ds.coords
    assert set(ds.coords['nb_faces'].values) == set([b.mesh.nb_faces for b in bodies])


def test_dataset_from_bemio():
    bemio = pytest.importorskip("bemio.io.wamit", reason="Bemio not installed, test skipped.")
    current_file_path = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(current_file_path, "Bemio_verification_cases", "sphere.out")
    bemio_data = bemio.read(out_file)

    new_dataset = cpt.assemble_dataset(bemio_data)
    assert (np.moveaxis(bemio_data.body[0].am.all, 2, 0) * bemio_data.body[0].rho == \
        new_dataset['added_mass'].values).all()
