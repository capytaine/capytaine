import numpy as np
import re
import os
from pathlib import Path
import xarray as xr
import capytaine as cpt
from capytaine.io.xarray import separate_complex_values
from capytaine.io.wamit import export_to_wamit
import random


def test_simulation_omega_zero_export():
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
    # Only keep data lines (skip headers/comments)
    data_lines = [
        line
        for line in lines
        if re.match(r"^-?\d+\.?\d*[eE]?[-+]?\d*\s+\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    assert data_lines, "No data lines found in the .1 file for omega=0."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 4
        ), f"For omega=0, each line in the .1 file must have 4 fields: {fields}"
    # Check the number of data lines: len(omega) * (6*6) = 1*36 = 36
    expected_lines = len(test_matrix["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."


def test_simulation_omega_inf_export():
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
    # Only keep data lines (skip headers/comments)
    data_lines = [
        line
        for line in lines
        if re.match(r"^-?\d+\.?\d*[eE]?[-+]?\d*\s+\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    assert data_lines, "No data lines found in the .1 file for omega=inf."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 4
        ), f"For omega=inf, each line in the .1 file must have 4 fields: {fields}"
    # Check the number of data lines: len(omega) * (6*6) = 1*36 = 36
    expected_lines = len(test_matrix["omega"]) * 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."


def test_simulation_omega_classic_export_hydrostatics():
    """
    Test a simulation with multiple classic omega values, export hydrostatics to WAMIT, and check the .hst file:
    - The file must exist
    - All data lines must have 4 fields
    - The number of data lines must be 1 (one line for hydrostatics)
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
    # Only keep data lines (skip headers/comments)
    data_lines = [
        line for line in lines if re.match(r"^\s*\d+\s+\d+\s+[-+eE0-9.]+$", line)
    ]
    print(data_lines)
    assert data_lines, "No data lines found in the .hst file for omega classic."
    for line in data_lines:
        fields = re.split(r"\s+", line)
        assert (
            len(fields) == 3
        ), f"For hydrostatics, each line in the .hst file must have 3 fields: {fields}"
    # Check the number of data lines: only one line for hydrostatics
    expected_lines = 36  # 6 dofs * 6 directions = 36
    assert (
        len(data_lines) == expected_lines
    ), f"Expected {expected_lines} data lines, found {len(data_lines)}."
