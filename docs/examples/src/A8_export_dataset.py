import numpy as np
from pathlib import Path
import xarray as xr

import capytaine as cpt
from capytaine.io.xarray import separate_complex_values
from capytaine.io.wamit import export_to_wamit

# --- Parameters ---
output_dir = Path("outputs")
output_dir.mkdir(exist_ok=True)

# --- Load mesh and set up body ---
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

# --- Export Capytaine NetCDF (.nc) ---
separate_complex_values(dataset).to_netcdf(
    output_dir / "test_boat200.nc",
    encoding={"radiating_dof": {"dtype": "U"}, "influenced_dof": {"dtype": "U"}},
)

# --- Export all WAMIT formats: .1, .3, .3fk, .3sc, .hst ---
export_to_wamit(
    dataset,
    problem_name=str(output_dir / "test_boat200"),
    exports=("1", "3", "3fk", "3sc", "hst"),
)

print("Export complete:")
print(f"- Capytaine NetCDF: {output_dir / 'test_boat200.nc'}")
print(f"- WAMIT .1 file: {output_dir / 'test_boat200.1'}")
print(f"- WAMIT .3 file: {output_dir / 'test_boat200.3'}")
print(f"- WAMIT .3fk file: {output_dir / 'test_boat200.3fk'}")
print(f"- WAMIT .3sc file: {output_dir / 'test_boat200.3sc'}")
print(f"- WAMIT .hst file: {output_dir / 'test_boat200.hst'}")
