import numpy as np
import xarray as xr
import capytaine as cpt
import unittest
import re
#from capytaine.io.xarray import separate_complex_values
from capytaine.io.export_wamit import export_to_wamit
from pathlib import Path

# =====================
# 1. Simulation
# =====================

# Load mesh and setup body
full_body = cpt.FloatingBody(
    mesh=cpt.load_mesh("docs/examples/src/boat_200.mar", file_format="nemoh"),
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
)
immersed_body = full_body.immersed_part()

# Setup simulation parameters
test_matrix = xr.Dataset({
    "omega": np.linspace(0.1, 1.0, 5),  # frequencies in rad/s.
    "wave_direction": np.linspace(0, np.pi, 3), # incident wave angles in radians.
    "radiating_dof": list(immersed_body.dofs), # degrees of freedom that radiate waves (Max 6)
    "water_depth": [np.inf], 
    "rho": [1025], # Water density in kg/m³
})

# Run simulation
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, immersed_body)

# Export to WAMIT format
export_dir = Path("my_examples")
export_dir.mkdir(exist_ok=True)
export_to_wamit(dataset, problem_name=str(export_dir / "boat_200"), exports=("1", "2", "3"))

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
        self.check_wamit_file("my_examples/boat_200.1", expected_fields=5, header_lines=1)

    def test_wamit_2(self):
        self.check_wamit_file("my_examples/boat_200.2", expected_fields=7, header_lines=1)

    def test_wamit_3(self):
        self.check_wamit_file("my_examples/boat_200.3", expected_fields=7, header_lines=1)


if __name__ == "__main__":
    print("\n✅ Running simulation and exporting to WAMIT...\n")
    unittest.main()
