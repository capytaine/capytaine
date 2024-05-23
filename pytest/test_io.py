import xarray as xr
import numpy as np
import capytaine as cpt
import pytest
import os

from capytaine.io.xarray import problems_from_dataset, separate_complex_values, merge_complex_values
from capytaine.io.legacy import export_hydrostatics


def test_incomplete_test_matrix():
    body = cpt.FloatingBody(cpt.mesh_sphere())
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
        cpt.FloatingBody(cpt.mesh_sphere(radius=1, resolution=(3, 3)), name="sphere_1"),
        cpt.FloatingBody(cpt.mesh_sphere(radius=2, resolution=(5, 3)), name="sphere_2"),
        cpt.FloatingBody(cpt.mesh_sphere(radius=3, resolution=(7, 3)), name="sphere_3"),
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


def test_legacy_export_hydrostatics():
    # create two cylinders & compute inertia & hydrostatic data
    cylinder1_mesh = cpt.meshes.predefined.mesh_vertical_cylinder(length=5.0, center = (1.0, 1.0, -0.5))
    cylinder1 = cpt.FloatingBody(cylinder1_mesh, center_of_mass=[1.0, 1.0, -2.0])
    cylinder1.add_all_rigid_body_dofs()
    cylinder1.compute_rigid_body_inertia()
    cylinder1.compute_hydrostatics()

    cylinder2_mesh = cpt.meshes.predefined.mesh_vertical_cylinder(length=5.0, center = (5.0, 1.0, -0.5))
    cylinder2 = cpt.FloatingBody(cylinder2_mesh, center_of_mass=[5.0, 1.0, -2.0])
    cylinder2.add_all_rigid_body_dofs()
    cylinder2.compute_rigid_body_inertia()
    cylinder2.compute_hydrostatics()

    current_file_path = os.path.dirname(os.path.abspath(__file__))
    testing_directory = os.path.join(current_file_path, 'io_legacy_cases')

    # export hydrostatics
    export_hydrostatics(os.path.join(testing_directory, 'single_body'), cylinder1)
    export_hydrostatics(os.path.join(testing_directory, 'single_body_list'), [cylinder1])
    export_hydrostatics(os.path.join(testing_directory, 'two_bodies_list'), [cylinder1, cylinder2])


    # Check single body Hydrostatics.dat
    with open(os.path.join(testing_directory, "single_body/Hydrostatics.dat"), "r") as f:
        single_body_Hydrostatics = f.read()
    with open(os.path.join(testing_directory, "reference_data/single_body/Hydrostatics.dat"), "r") as f:
        single_body_Hydrostatics_ref = f.read()
    assert single_body_Hydrostatics == single_body_Hydrostatics_ref


    # Check single body (list) Hydrostatics.dat
    with open(os.path.join(testing_directory, "single_body_list/Hydrostatics.dat"), "r") as f:
        single_body_list_Hydrostatics = f.read()
    with open(os.path.join(testing_directory, "reference_data/single_body_list/Hydrostatics.dat"), "r") as f:
        single_body_list_Hydrostatics_ref = f.read()
    assert single_body_list_Hydrostatics == single_body_list_Hydrostatics_ref


    # Check two bodies (list) Hydrostatics_0.dat
    with open(os.path.join(testing_directory, "two_bodies_list/Hydrostatics_0.dat"), "r") as f:
        two_bodies_Hydrostatics_0 = f.read()
    with open(os.path.join(testing_directory, "reference_data/two_bodies_list/Hydrostatics_0.dat"), "r") as f:
        two_bodies_Hydrostatics_0_ref = f.read()
    assert two_bodies_Hydrostatics_0 == two_bodies_Hydrostatics_0_ref

    # Check two bodies (list) Hydrostatics_1.dat
    with open(os.path.join(testing_directory, "two_bodies_list/Hydrostatics_1.dat"), "r") as f:
        two_bodies_Hydrostatics_1 = f.read()
    with open(os.path.join(testing_directory, "reference_data/two_bodies_list/Hydrostatics_1.dat"), "r") as f:
        two_bodies_Hydrostatics_1_ref = f.read()
    assert two_bodies_Hydrostatics_1 == two_bodies_Hydrostatics_1_ref


    # Check single body KH.dat
    single_body_KH = np.loadtxt(os.path.join(testing_directory, "single_body/KH.dat"))
    single_body_KH_ref = np.loadtxt(os.path.join(testing_directory, "reference_data/single_body/KH.dat"))
    np.testing.assert_allclose(single_body_KH, single_body_KH_ref, atol=1e-6)

    # Check single body (list) KH.dat
    single_body_list_KH = np.loadtxt(os.path.join(testing_directory, "single_body_list/KH.dat"))
    single_body_list_KH_ref = np.loadtxt(os.path.join(testing_directory, "reference_data/single_body_list/KH.dat"))
    np.testing.assert_allclose(single_body_list_KH, single_body_list_KH_ref, atol=1e-6)

    # Check two bodies (list) KH_0.dat
    two_bodies_KH_0 = np.loadtxt(os.path.join(testing_directory, "two_bodies_list/KH_0.dat"))
    two_bodies_KH_0_ref = np.loadtxt(os.path.join(testing_directory, "reference_data/two_bodies_list/KH_0.dat"))
    np.testing.assert_allclose(two_bodies_KH_0, two_bodies_KH_0_ref, atol=1e-6)

    # Check two bodies (list) KH_1.dat
    two_bodies_KH_1 = np.loadtxt(os.path.join(testing_directory, "two_bodies_list/KH_1.dat"))
    two_bodies_KH_1_ref = np.loadtxt(os.path.join(testing_directory, "reference_data/two_bodies_list/KH_1.dat"))
    np.testing.assert_allclose(two_bodies_KH_1, two_bodies_KH_1_ref, atol=1e-6)
