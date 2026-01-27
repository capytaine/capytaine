import numpy as np
import xarray as xr
import capytaine as cpt

cpt.set_logging("INFO")

# Generating a geometric mesh.
mesh = cpt.mesh_horizontal_cylinder(
    length=10.0,
    radius=1.0,
    center=(0, 0, 0,),
    resolution=(10, 20, 30,),
    name="cylinder_mesh",
)
# Consider also using `cpt.load_mesh` to open an existing file.

# Define a rigid body using this mesh
body = cpt.FloatingBody(
    mesh=mesh,
    lid_mesh=mesh.generate_lid(),  # Generate a lid for irregular frequency removal
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
    name="floating cylinder",
)

# Set up the problems: we will solve a radiation problem for each
# degree of freedom and a diffraction problem for each wave direction,
# both for each frequency in the frequency range.
test_matrix = xr.Dataset(
    {
        "omega": np.linspace(0.1, 3.0, 20),
        "wave_direction": [0.0, np.pi/2],
        "radiating_dof": list(body.dofs),
    }
)
# Water density, gravity and water depth have not been specified.
# Default values are used.

# Solve all radiation problems
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, body.immersed_part())


# Export data in various formats
from pathlib import Path
output_dir = Path("A1_example_outputs")
output_dir.mkdir(exist_ok=True)
cpt.export_dataset(output_dir / "test_boat200.nc", dataset)  # Format is deduced from filename suffix.
cpt.export_dataset(output_dir / "wamit_data", dataset, format="wamit")
cpt.export_dataset(output_dir / "nemoh_data", dataset, format="nemoh")


# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
for dof in body.dofs:
    plt.plot(
        dataset.coords["omega"],
        dataset["added_mass"].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker="o",
    )
plt.xlabel("omega")
plt.ylabel("added mass")
plt.legend()
plt.tight_layout()
plt.show()
