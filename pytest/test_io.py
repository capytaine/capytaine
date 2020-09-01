
import xarray as xr
import numpy as np
import capytaine as cpt

from capytaine.io.xarray import separate_complex_values, merge_complex_values


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
