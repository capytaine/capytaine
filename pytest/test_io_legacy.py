import os

from functools import lru_cache
import numpy as np
import pytest

import capytaine as cpt
from capytaine.io.legacy import export_hydrostatics


@lru_cache
def setup_two_bodies():
    cylinder1_mesh = cpt.mesh_vertical_cylinder(length=5.0, center=(1.0, 1.0, -0.5))
    cylinder1 = cpt.FloatingBody(
        cylinder1_mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(1.0, 1.0, -2.0)),
        center_of_mass=[1.0, 1.0, -2.0],
        name="Cylinder1"
    )

    cylinder2_mesh = cpt.mesh_vertical_cylinder(length=5.0, center=(5.0, 1.0, -0.5))
    cylinder2 = cpt.FloatingBody(
        cylinder2_mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 1.0, -2.0)),
        center_of_mass=[5.0, 1.0, -2.0],
        name="Cylinder2"
    )
    return cylinder1, cylinder2

def setup_two_bodies_with_precomputed_hydrostatics():
    cylinder1, cylinder2 = setup_two_bodies()
    cylinder1.hydrostatic_stiffness = cylinder1.immersed_part().compute_hydrostatic_stiffness()
    cylinder2.hydrostatic_stiffness = cylinder2.immersed_part().compute_hydrostatic_stiffness()
    return cylinder1, cylinder2

@lru_cache
def load_reference_data_single_body():
    current_file_path = os.path.dirname(os.path.abspath(__file__))
    data_directory = os.path.join(current_file_path, 'io_legacy_cases')
    with open(os.path.join(data_directory, "reference_data/single_body/Hydrostatics.dat"), "r") as f:
       hydrostatics = f.read()
    KH = np.loadtxt(os.path.join(data_directory, "reference_data/single_body/KH.dat"))
    return {"Hydrostatics.dat": hydrostatics, "KH.dat": KH}

def check_single_body_hydrostatics_content(file_path):
    ref_data = load_reference_data_single_body()

    with open(os.path.join(file_path, "Hydrostatics.dat"), "r") as f:
        single_body_Hydrostatics = f.read()
    assert single_body_Hydrostatics == ref_data["Hydrostatics.dat"]

    single_body_KH = np.loadtxt(os.path.join(file_path, "KH.dat"))
    np.testing.assert_allclose(single_body_KH, ref_data["KH.dat"], atol=1e-6)

@lru_cache
def load_reference_data_two_bodies():
    current_file_path = os.path.dirname(os.path.abspath(__file__))
    data_directory = os.path.join(current_file_path, 'io_legacy_cases')
    with open(os.path.join(data_directory, "reference_data/two_bodies/Hydrostatics_0.dat"), "r") as f:
        hydrostatics_0 = f.read()
    with open(os.path.join(data_directory, "reference_data/two_bodies/Hydrostatics_1.dat"), "r") as f:
        hydrostatics_1 = f.read()
    KH_0 = np.loadtxt(os.path.join(data_directory, "reference_data/two_bodies/KH_0.dat"))
    KH_1 = np.loadtxt(os.path.join(data_directory, "reference_data/two_bodies/KH_1.dat"))
    return {"Hydrostatics_0.dat": hydrostatics_0, "Hydrostatics_1.dat": hydrostatics_1, "KH_0.dat": KH_0, "KH_1.dat": KH_1}

def check_two_bodies_hydrostatics_content(file_path):
    ref_data = load_reference_data_two_bodies()

    with open(os.path.join(file_path, "Hydrostatics_0.dat"), "r") as f:
        two_bodies_Hydrostatics_0 = f.read()
    assert two_bodies_Hydrostatics_0 == ref_data["Hydrostatics_0.dat"]
    with open(os.path.join(file_path, "Hydrostatics_1.dat"), "r") as f:
        two_bodies_Hydrostatics_1 = f.read()
    assert two_bodies_Hydrostatics_1 == ref_data["Hydrostatics_1.dat"]

    two_bodies_KH_0 = np.loadtxt(os.path.join(file_path, "KH_0.dat"))
    np.testing.assert_allclose(two_bodies_KH_0, ref_data["KH_0.dat"], atol=1e-6)
    two_bodies_KH_1 = np.loadtxt(os.path.join(file_path, "KH_1.dat"))
    np.testing.assert_allclose(two_bodies_KH_1, ref_data["KH_1.dat"], atol=1e-6)


def test_legacy_export_hydrostatics_single_body(tmp_path):
    cylinder1, _ = setup_two_bodies_with_precomputed_hydrostatics()
    export_hydrostatics(tmp_path, cylinder1)
    check_single_body_hydrostatics_content(tmp_path)

def test_legacy_export_hydrostatics_single_body_as_list(tmp_path):
    cylinder1, _ = setup_two_bodies_with_precomputed_hydrostatics()
    export_hydrostatics(tmp_path, [cylinder1])
    check_single_body_hydrostatics_content(tmp_path)

def test_legacy_export_hydrostatics_two_bodies_as_list(tmp_path):
    cylinder1, cylinder2 = setup_two_bodies_with_precomputed_hydrostatics()
    export_hydrostatics(tmp_path, [cylinder1, cylinder2])
    check_two_bodies_hydrostatics_content(tmp_path)

def test_legacy_export_hydrostatics_two_bodies_as_multibody(tmp_path):
    cylinder1, cylinder2 = setup_two_bodies_with_precomputed_hydrostatics()
    multibody = cpt.Multibody([cylinder1, cylinder2])
    export_hydrostatics(tmp_path, multibody)
    check_two_bodies_hydrostatics_content(tmp_path)

def test_legacy_export_hydrostatics_two_bodies_computing_hydrostatics_on_multibody(tmp_path):
    # Unlike the previous example, the hydrostatic stiffness is not precomputed on the individual bodies, but directly on the multibody.
    cylinder1, cylinder2 = setup_two_bodies()
    multibody = cpt.Multibody([cylinder1, cylinder2])
    multibody.hydrostatic_stiffness = multibody.immersed_part().compute_hydrostatic_stiffness()
    with pytest.raises(ValueError):
        export_hydrostatics(tmp_path, multibody)
