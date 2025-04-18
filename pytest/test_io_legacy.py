import os

import numpy as np
import pytest

import capytaine as cpt
from capytaine.io.legacy import export_hydrostatics


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
