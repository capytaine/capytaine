import numpy as np
import xarray as xr
import capytaine as cpt

body_1 = cpt.FloatingBody(
            mesh=cpt.load_mesh("./boat_200.mar", file_format="nemoh"),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        )
# If you have several rigid bodies, copy the code above to define "body_2", "body_3", etc.

all_bodies = body_1  # Replace "body_1" by "body_1 + body_2 + body_3" for multibody problem.

all_bodies = all_bodies.immersed_part()  # if your mesh has panels above the free surface, remove them

# Set up paramaters
test_matrix = xr.Dataset({
        "omega": np.linspace(0.1, 2.0, 20),  # Can also specify "period", "wavelength" or "wavenumber"
        "wave_direction": np.linspace(0, np.pi, 3),
        "radiating_dof": list(all_bodies.dofs),
        "water_depth": [np.inf],
        "rho": [1025],
        })

# Do the resolution
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, all_bodies)

# Export to NetCDF file
from capytaine.io.xarray import separate_complex_values
separate_complex_values(dataset).to_netcdf("dataset.nc",
                                           encoding={'radiating_dof': {'dtype': 'U'},
                                                     'influenced_dof': {'dtype': 'U'}})
