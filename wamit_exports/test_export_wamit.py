import numpy as np
import unittest
import re
import os
from pathlib import Path
import xarray as xr
import capytaine as cpt
from capytaine.io.xarray import separate_complex_values
from capytaine.io.wamit import export_to_wamit

# ================================
# 1. Simulation
# ================================

# Load mesh and setup body
mesh = cpt.load_mesh("docs/examples/src/boat_200.mar", file_format="nemoh")
dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
full_body = cpt.FloatingBody(mesh, dofs)
immersed_body = full_body.immersed_part()

# Setup simulation parameters
test_matrix = xr.Dataset({
    # Incident wave frequencies in rad/s
    "omega": [0] + np.linspace(0.1, 1.0, 5).tolist() + [np.inf],
    # Incident wave angles in rad
    "wave_direction": np.linspace(0, np.pi, 3),
    # Degrees of freedom that radiate waves (['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])
    "radiating_dof": list(immersed_body.dofs),
    # Water depth (float or np.inf)
    "water_depth": [np.inf], 
    # Water density in kg/m³
    "rho": [1025],
})


# Run simulation using parameters matrix
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, immersed_body)

# Export simulation results to ...
export_dir = Path("wamit_exports")
export_dir.mkdir(exist_ok=True)
# ... NetCDF format
separate_complex_values(dataset).to_netcdf(
    export_dir/"boat_200.nc",
    encoding={
        'radiating_dof': {'dtype': 'U'},
        'influenced_dof': {'dtype': 'U'}
    }
)
# ... WAMIT format
export_to_wamit(dataset, problem_name=str(export_dir / "boat_200"), 
                exports=("1", "3", "3fk", "3sc"))


# =====================
# 2. Unit Tests
# =====================

class TestWAMITFileFormat(unittest.TestCase):
    def check_wamit_file(self, filepath, expected_fields, header_lines=1):
        with open(filepath, "r") as f:
            lines = f.readlines()[header_lines:]

        for line in lines:
            if not line.strip():
                continue
            fields = re.split(r"\s+", line.strip())
            if len(fields) in (1, 2):  # index or DOF header lines
                continue
            self.assertEqual(
                len(fields), expected_fields,
                msg=f"Line has {len(fields)} fields, expected {expected_fields}: {line.strip()}"
            )
            for val in fields:
                try:
                    float(val)
                except ValueError:
                    self.fail(f"Non-numeric value found: {val} in line: {line.strip()}")

    def test_wamit_1(self):
        self.check_wamit_file("wamit_exports/boat_200.1", expected_fields=5, header_lines=1)

    def test_wamit_3(self):
        self.check_wamit_file("wamit_exports/boat_200.3", expected_fields=7, header_lines=1)

    def test_wamit_3fk(self):
        self.check_wamit_file("wamit_exports/boat_200.3fk", expected_fields=7, header_lines=1)

    def test_wamit_3sc(self):
        self.check_wamit_file("wamit_exports/boat_200.3sc", expected_fields=7, header_lines=1)

class TestWAMITPhysicalConsistency(unittest.TestCase):
    def setUp(self):
        self.base_filename = "wamit_exports/boat_200"
        self.nc_path = self.base_filename + ".nc"
        self.wamit_1_path = self.base_filename + ".1"
        self.wamit_3_path = self.base_filename + ".3"
        self.wamit_3fk_path = self.base_filename + ".3fk"
        self.wamit_3sc_path = self.base_filename + ".3sc"

        # Load dataset
        self.dataset = xr.open_dataset(self.nc_path)

        # Read all exported WAMIT text files
        self.data_1 = self._read_file(self.wamit_1_path)
        self.data_3 = self._read_file(self.wamit_3_path, skip_comments=True)
        self.data_3fk = self._read_file(self.wamit_3fk_path, skip_comments=True)
        self.data_3sc = self._read_file(self.wamit_3sc_path, skip_comments=True)

    def _read_file(self, path, skip_comments=False):
        with open(path, "r") as f:
            return [
                line.strip() for line in f
                if line.strip() and (not skip_comments or not line.startswith("#"))
            ]

    def _extract_periods(self, lines, expected_fields):
        periods = set()
        for line in lines:
            fields = re.split(r'\s+', line.strip())
            if len(fields) == expected_fields:
                try:
                    periods.add(round(float(fields[0]), 6))
                except ValueError:
                    continue
        return periods

    def test_exported_files_exist(self):
        for path in [
            self.nc_path,
            self.wamit_1_path,
            self.wamit_3_path,
            self.wamit_3fk_path,
            self.wamit_3sc_path,
        ]:
            self.assertTrue(os.path.isfile(path), f"Missing file: {path}")

    def test_frequencies_match_in_files(self):
        expected_periods = np.round(2 * np.pi / self.dataset["omega"].values, 6)

        periods_1 = self._extract_periods(self.data_1, expected_fields=5)
        periods_3 = self._extract_periods(self.data_3, expected_fields=7)
        periods_3fk = self._extract_periods(self.data_3fk, expected_fields=7)
        periods_3sc = self._extract_periods(self.data_3sc, expected_fields=7)

        for per in expected_periods:
            self.assertTrue(
                any(np.isclose(per, p, atol=1e-5) for p in periods_1),
                f"Period {per:.6f} missing in .1 file"
            )
            self.assertTrue(
                any(np.isclose(per, p, atol=1e-5) for p in periods_3),
                f"Period {per:.6f} missing in .3 file"
            )
            self.assertTrue(
                any(np.isclose(per, p, atol=1e-5) for p in periods_3fk),
                f"Period {per:.6f} missing in .3fk file"
            )
            self.assertTrue(
                any(np.isclose(per, p, atol=1e-5) for p in periods_3sc),
                f"Period {per:.6f} missing in .3sc file"
            )
            
if __name__ == "__main__":
    print("\n✅ Running simulation and exporting to WAMIT...\n")
    unittest.main()
