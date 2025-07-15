import numpy as np
from pathlib import Path
import xarray as xr

import capytaine as cpt

# --- Parameters ---
output_dir = Path("outputs")
output_dir.mkdir(exist_ok=True)

# --- Create mesh and set up body ---
mesh = cpt.mesh_sphere()
dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
full_body = cpt.FloatingBody(mesh, dofs)
full_body.center_of_mass = np.copy(full_body.mesh.center_of_mass_of_nodes)
immersed_body = full_body.immersed_part()
immersed_body.compute_hydrostatics()

# --- Define parameter grid ---
omegas = np.concatenate(([0], np.linspace(0.1, 1.0, 10), [np.inf]))  # 0, 5 values, inf
wave_directions = np.linspace(0, np.pi, 3)  # 3 directions: 0, pi/2, pi

test_matrix = xr.Dataset(
    {
        "omega": omegas,
        "wave_direction": wave_directions,
        "radiating_dof": list(immersed_body.dofs),
        "water_depth": [np.inf],
        "rho": [1025],
    }
)

# --- Run simulation ---
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, immersed_body)


cpt.export_dataset(output_dir / "test_boat200.nc", dataset)  # Format is deduced from filename suffix.
cpt.export_dataset(output_dir / "wamit_data", dataset, format="wamit")
(output_dir / 'nemoh_data').mkdir(exist_ok=True)
cpt.export_dataset(output_dir / "nemoh_data", dataset, format="nemoh")

print("Exported files:")
for file in output_dir.glob("**/*"):
    print(file)
