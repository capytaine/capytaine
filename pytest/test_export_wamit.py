import os
import numpy as np
import re
from pathlib import Path
import xarray as xr
import capytaine as cpt
from capytaine.io.wamit import export_to_wamit
import random
import pandas as pd


def test_wamit_export_and_import_force():
    """
    Run a capytaine simulation for a unique known frequency, export to WAMIT (.1 and .3 files),
    then reload and check that the frequency matches.
    """
    mar_file = (
        Path(__file__).parent.parent / "docs" / "examples" / "src" / "boat_200.mar"
    )
    mesh = cpt.load_mesh(str(mar_file), file_format="nemoh")
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()
    immersed_body.compute_hydrostatics()
    known_omega = 0.75
    known_period = 2 * np.pi / known_omega
    test_matrix = xr.Dataset(
        {
            "omega": [known_omega],
            "wave_direction": [0.0],
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )
    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(test_matrix, immersed_body)
    export_dir = Path(__file__).parent.parent / "pytest"
    problem_name = str(export_dir / "boat_200_unique")
    export_to_wamit(dataset, problem_name=problem_name, exports=("1", "3"))
    wamit_1_path = export_dir / "boat_200_unique.1"
    wamit_3_path = export_dir / "boat_200_unique.3"
    assert wamit_1_path.exists(), f"File {wamit_1_path} was not generated."
    assert wamit_3_path.exists(), f"File {wamit_3_path} was not generated."
    # Now import the .1 file and check the frequency
    imported_df = pd.read_csv(wamit_1_path, sep=r"\s+", header=None)
    # The first column is PERIOD
    imported_period = float(imported_df.iloc[0, 0])
    assert np.isclose(
        imported_period, known_period
    ), f"Imported period {imported_period} does not match expected {known_period}"
    os.remove(wamit_1_path)

    # Now import the .3 file and check the frequency
    imported_df = pd.read_csv(wamit_3_path, sep=r"\s+", header=None)
    # The first column is omega
    imported_period = float(imported_df.iloc[0, 0])
    assert np.isclose(
        imported_period, known_period
    ), f"Imported period {imported_period} does not match expected {known_period}"
    os.remove(wamit_3_path)


def test_wamit_export_omega_zero():
    """
    Test a simulation with omega=0, export to WAMIT, and check the .1 file:
    - The file must exist
    - All data lines must have 4 fields
    - The number of data lines must be len(omega) * 36 (here: 36)
    """
    mar_file = (
        Path(__file__).parent.parent / "docs" / "examples" / "src" / "boat_200.mar"
    )
    mesh = cpt.load_mesh(str(mar_file), file_format="nemoh")
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()
    immersed_body.compute_hydrostatics()
    test_matrix = xr.Dataset(
        {
            "omega": [0],
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )
    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(test_matrix, immersed_body)
    export_dir = Path(__file__).parent.parent / "pytest"
    export_to_wamit(
        dataset, problem_name=str(export_dir / "boat_200_omega0"), exports=("1",)
    )
    wamit_1_path = export_dir / "boat_200_omega0.1"
    assert wamit_1_path.exists(), f"File {wamit_1_path} was not generated."
    with open(wamit_1_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    data_lines = [
        line
        for line in lines
        if re.match(r"^\s*[-+eE0-9.]+\s+\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    assert data_lines, "No data lines found in the .1 file for omega=0."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 4
        ), f"For omega=0, each line in the .1 file must have 4 fields: {fields}"
    expected_lines = len(test_matrix["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."
    os.remove(wamit_1_path)  # Clean up after test


def test_wamit_export_omega_inf():
    """
    Test a simulation with omega=inf, export to WAMIT, and check the .1 file:
    - The file must exist
    - All data lines must have 4 fields
    - The number of data lines must be len(omega) * 36 (here: 36)
    """
    mar_file = (
        Path(__file__).parent.parent / "docs" / "examples" / "src" / "boat_200.mar"
    )
    mesh = cpt.load_mesh(str(mar_file), file_format="nemoh")
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()
    immersed_body.compute_hydrostatics()
    test_matrix = xr.Dataset(
        {
            "omega": [np.inf],
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )
    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(test_matrix, immersed_body)
    export_dir = Path(__file__).parent
    export_to_wamit(
        dataset, problem_name=str(export_dir / "boat_200_omegaINF"), exports=("1",)
    )
    wamit_1_path = export_dir / "boat_200_omegaINF.1"
    assert wamit_1_path.exists(), f"File {wamit_1_path} was not generated."
    with open(wamit_1_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    data_lines = [
        line
        for line in lines
        if re.match(r"^\s*[-+eE0-9.]+\s+\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    assert data_lines, "No data lines found in the .1 file for omega=inf."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 4
        ), f"For omega=inf, each line in the .1 file must have 4 fields: {fields}"
    expected_lines = len(test_matrix["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."
    os.remove(wamit_1_path)  # Clean up after test


def test_wamit_export_hydrostatics():
    """
    Test a simulation with multiple classic omega values, export hydrostatics to WAMIT, and check the .hst file:
    - The file must exist
    - All data lines must have 3 fields
    - The number of data lines must be 36 (6x6 matrix, possibly sparse)
    """
    mar_file = (
        Path(__file__).parent.parent / "docs" / "examples" / "src" / "boat_200.mar"
    )
    mesh = cpt.load_mesh(str(mar_file), file_format="nemoh")
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()
    immersed_body.compute_hydrostatics()
    test_matrix = xr.Dataset(
        {
            "omega": np.linspace(0.1, 1.0, 5),
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )
    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(test_matrix, immersed_body)
    export_dir = Path(__file__).parent.parent / "pytest"
    export_to_wamit(
        dataset, problem_name=str(export_dir / "boat_200_omega"), exports=("hst",)
    )
    wamit_hst_path = export_dir / "boat_200_omega.hst"
    assert wamit_hst_path.exists(), f"File {wamit_hst_path} was not generated."
    with open(wamit_hst_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    data_lines = [
        line for line in lines if re.match(r"^\s*\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    assert data_lines, "No data lines found in the .hst file for omega classic."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 3
        ), f"For hydrostatics, each line in the .hst file must have 3 fields: {fields}"
    expected_lines = 36  # 6 dofs * 6 directions = 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."
    os.remove(wamit_hst_path)  # Clean up after test
