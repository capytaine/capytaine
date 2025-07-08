import re
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import pytest

import capytaine as cpt
from capytaine.io.wamit import export_to_wamit


@pytest.fixture
def full_dataset():
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
            "omega": [0.0, 1.0, np.inf],
            "wave_direction": [0.0],
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )
    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(test_matrix, immersed_body)
    return dataset


def test_export_wamit_and_import_force(full_dataset, tmpdir):
    """Export to WAMIT (.1 and .3 files), then reload and check that the frequency matches."""
    dataset = full_dataset.sel(omega=[1.0])
    problem_name = str(tmpdir / "boat_200_unique")
    export_to_wamit(dataset, problem_name=problem_name, exports=("1", "3"))
    wamit_1_path = tmpdir / "boat_200_unique.1"
    wamit_3_path = tmpdir / "boat_200_unique.3"
    assert wamit_1_path.exists(), f"File {wamit_1_path} was not generated."
    assert wamit_3_path.exists(), f"File {wamit_3_path} was not generated."

    # Now import the .1 file and check the frequency
    imported_df = pd.read_csv(wamit_1_path, sep=r"\s+", header=None)
    # The first column is PERIOD
    imported_period = float(imported_df.iloc[0, 0])
    ref_period = float(dataset.period)
    assert np.isclose(
        imported_period, ref_period
    ), f"Imported period {imported_period} does not match expected {ref_period}"

    # Now import the .3 file and check the frequency
    imported_df = pd.read_csv(wamit_3_path, sep=r"\s+", header=None)
    # The first column is omega
    imported_period = float(imported_df.iloc[0, 0])
    assert np.isclose(
        imported_period, ref_period
    ), f"Imported period {imported_period} does not match expected {ref_period}"


def test_export_wamit_omega_zero(full_dataset, tmpdir):
    """
    Test a simulation with omega=0, export to WAMIT, and check the .1 file:
    - The file must exist
    - All data lines must have 4 fields
    - The number of data lines must be len(omega) * 36 (here: 36)
    """
    dataset = full_dataset.sel(omega=[0.0])
    export_to_wamit(
        dataset, problem_name=str(tmpdir / "boat_200_omega0"), exports=("1",)
    )
    wamit_1_path = tmpdir / "boat_200_omega0.1"
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
    expected_lines = len(dataset["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."


def test_export_wamit_omega_inf(full_dataset, tmpdir):
    """
    Test a simulation with omega=inf, export to WAMIT, and check the .1 file:
    - The file must exist
    - All data lines must have 4 fields
    - The number of data lines must be len(omega) * 36 (here: 36)
    """
    dataset = full_dataset.sel(omega=[np.inf])
    export_to_wamit(
        dataset, problem_name=str(tmpdir / "boat_200_omegaINF"), exports=("1",)
    )
    wamit_1_path = tmpdir / "boat_200_omegaINF.1"
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
    expected_lines = len(dataset["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."


def test_export_wamit_hydrostatics(full_dataset, tmpdir):
    """
    Test a simulation with multiple classic omega values, export hydrostatics to WAMIT, and check the .hst file:
    - The file must exist
    - All data lines must have 3 fields
    - The number of data lines must be 36 (6x6 matrix, possibly sparse)
    """
    export_to_wamit(
        full_dataset, problem_name=str(tmpdir / "boat_200_omega"), exports=("hst",)
    )
    wamit_hst_path = tmpdir / "boat_200_omega.hst"
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


def test_warning_for_wrong_format(full_dataset, tmpdir, caplog):
    """Test the warning when trying to export excitation force but they are not in the dataset"""
    radiation_only_ds = full_dataset[["added_mass", "radiation_damping"]]
    problem_name = str(tmpdir / 'radiation_data')
    with caplog.at_level(logging.WARNING):
        export_to_wamit(
                radiation_only_ds, problem_name=problem_name, exports=["3"]
                )
    assert "Missing field 'excitation_force'" in caplog.text
    assert not (tmpdir / 'radiation_data.3').exists()
