import numpy as np
import xarray as xr
import capytaine as cpt

mesh_1 = cpt.load_mesh("./boat_200.mar", file_format="nemoh")
body_1 = cpt.FloatingBody(
            mesh=mesh_1,
            lid_mesh=mesh_1.generate_lid(faces_max_radius=1.0),  # Adjust to the desired mesh resolution or remove lid entirely
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0)  # Optional, only for hydrostatics
        )
body_1.inertia_matrix = body_1.compute_rigid_body_inertia(rho=1025)
body_1.hydrostatic_stiffness = body_1.immersed_part().compute_hydrostatic_stiffness(rho=1025)

# If you have several rigid bodies, copy the code above to define "body_2", "body_3", etc.

all_bodies = body_1  # Replace "body_1" by "body_1 + body_2 + body_3" for multibody problem.
all_bodies = all_bodies.immersed_part()  # if the mesh has panels above the free surface, this command removes them

# Set up parameters
test_matrix = xr.Dataset({
        "omega": np.linspace(0.1, 2.0, 20),  # Can also specify "freq" (in Hz), "period", "wavelength" or "wavenumber"
        "wave_direction": np.linspace(0, np.pi, 3),
        "radiating_dof": list(all_bodies.dofs),
        "water_depth": [np.inf],
        "rho": [1025],
        })

# Do the resolution
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, all_bodies)

# Export to netcdf file
cpt.export_dataset("capy_dataset.nc", dataset)
