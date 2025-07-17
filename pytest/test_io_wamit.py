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
    mesh = cpt.mesh_sphere(resolution=(4, 4))
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


@pytest.fixture
def dataset_with_multiple_rho():
    mesh = cpt.mesh_sphere(resolution=(4, 4))
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()

    # Define a test matrix with multiple rho
    test_matrix = xr.Dataset(
        {
            "omega": ("omega", [1.0]),
            "wave_direction": ("wave_direction", [0.0]),
            "radiating_dof": ("radiating_dof", list(immersed_body.dofs)),
            "water_depth": ("water_depth", [np.inf]),
            "rho": ("rho", [1020.0, 1025.0]),  # ➜ should lead to an error
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


def test_export_wamit_hydrostatic(full_dataset, tmpdir):
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
    problem_name = str(tmpdir / "radiation_data")
    with caplog.at_level(logging.WARNING):
        export_to_wamit(radiation_only_ds, problem_name=problem_name, exports=["3"])
    assert "Missing field 'excitation_force'" in caplog.text
    assert not (tmpdir / "radiation_data.3").exists()


def compare_wamit_files(file1, file2, atol=1e-14):
    """Compare two WAMIT output files by parsing numerical values.

    Args:
        file1 (str): Path to the first WAMIT file.
        file2 (str): Path to the second WAMIT file.
        atol (float): Absolute tolerance for numerical comparison.

    Raises:
        AssertionError: If line counts differ or any value mismatches beyond tolerance.
    """
    with open(file1, "r", encoding="utf-8") as f:
        lines1 = [line.strip() for line in f if line.strip()]
    with open(file2, "r", encoding="utf-8") as f:
        lines2 = [line.strip() for line in f if line.strip()]

    assert len(lines1) == len(
        lines2
    ), f"Line count mismatch: {len(lines1)} vs {len(lines2)}"

    for i, (l1, l2) in enumerate(zip(lines1, lines2)):
        try:
            vals1 = [float(s) for s in re.split(r"\s+", l1)]
            vals2 = [float(s) for s in re.split(r"\s+", l2)]
        except ValueError:
            raise AssertionError(f"Line {i} contains non-numeric data:\n{l1}\n{l2}")

        assert len(vals1) == len(
            vals2
        ), f"Line {i} has different number of columns: {vals1} vs {vals2}"

        def condition(vals1, vals2, atol):
            return not np.allclose(vals1, vals2, atol=atol)

        if condition(vals1, vals2, atol):
            diffs = [
                f"{v1} ≠ {v2}"
                for v1, v2 in zip(vals1, vals2)
                if condition([v1], [v2], atol)
            ]
            raise AssertionError(
                f"Line {i} values differ beyond tolerance: {l1} // {l2} -- Differences: {', '.join(diffs)}"
            )


rng = np.random.default_rng(seed=42)
omega_values = rng.uniform(0.5, 5.0, size=3)


def test_export_wamit_fails_with_multiple_rho(dataset_with_multiple_rho, tmpdir):
    """Test that export fails when multiple rho values are present."""
    problem_name = tmpdir / "boat_invalid_rho"

    with pytest.raises(ValueError, match="only one value.*rho"):
        export_to_wamit(dataset_with_multiple_rho, problem_name=str(problem_name))


@pytest.mark.parametrize("omega_val", omega_values)
@pytest.mark.parametrize("export_type", ["hst", "1", "3"])
def test_export_wamit_frequency_axis_representations(export_type, omega_val, tmpdir):
    """
    Run Capytaine simulation with a single omega, then:
    - create two datasets: one with omega as dimension, one with period
    - export both
    - compare the resulting .1 and .3 files
    """
    # Load mesh and run BEM simulation for 1 omega
    mesh = cpt.mesh_sphere(resolution=(4, 4))
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
    full_body = cpt.FloatingBody(mesh, dofs)
    full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
    immersed_body = full_body.immersed_part()
    if export_type == "hst":
        immersed_body.compute_hydrostatics()

    test_matrix_omega = xr.Dataset(
        {
            "omega": np.asarray(omega_val),
            "wave_direction": [0.0],
            "radiating_dof": list(immersed_body.dofs),
            "water_depth": [np.inf],
            "rho": [1025],
        }
    )

    solver = cpt.BEMSolver()
    ds_omega = solver.fill_dataset(test_matrix_omega, immersed_body)

    # Create period-based version
    ds_period = ds_omega.swap_dims({"omega": "period"})
    ds_wavenumber = ds_period.swap_dims({"period": "wavenumber"})
    ds_wavelength = ds_wavenumber.swap_dims({"wavenumber": "wavelength"})

    # Export all
    tmpdir = Path(tmpdir)
    omega_name = tmpdir / f"omega_export_{export_type}"
    period_name = tmpdir / f"period_export_{export_type}"
    wavenumber_name = tmpdir / f"wavenumber_export_{export_type}"
    wavelength_name = tmpdir / f"wavelength_export_{export_type}"

    export_to_wamit(ds_omega, problem_name=str(omega_name), exports=[export_type])
    export_to_wamit(ds_period, problem_name=str(period_name), exports=[export_type])
    export_to_wamit(
        ds_wavenumber, problem_name=str(wavenumber_name), exports=[export_type]
    )
    export_to_wamit(
        ds_wavelength, problem_name=str(wavelength_name), exports=[export_type]
    )

    f_omega = omega_name.with_name(f"{omega_name.name}.{export_type}")
    f_period = period_name.with_name(f"{period_name.name}.{export_type}")
    f_wavenumber = wavenumber_name.with_name(f"{wavenumber_name.name}.{export_type}")
    f_wavelength = wavelength_name.with_name(f"{wavelength_name.name}.{export_type}")

    assert f_omega.exists(), f"File {f_omega} not generated"
    assert f_period.exists(), f"File {f_period} not generated"
    assert f_wavenumber.exists(), f"File {f_wavenumber} not generated"
    assert f_wavelength.exists(), f"File {f_wavelength} not generated"

    compare_wamit_files(f_omega, f_period)
    compare_wamit_files(f_omega, f_wavenumber)
    compare_wamit_files(f_omega, f_wavelength)
