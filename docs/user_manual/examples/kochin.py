import numpy as np
from numpy import pi
import xarray as xr
import matplotlib.pyplot as plt

import capytaine as cpt

cpt.set_logging('INFO')

# SET UP THE MESH
buoy = cpt.FloatingBody(
        mesh=cpt.mesh_vertical_cylinder(
            radius=1.0, length=3.0, center=(0, 0, -1.0),
            resolution=(15, 15, 3)
            ))
buoy.add_translation_dof(name="Surge")

# SOLVE THE BEM PROBLEMS AND COMPUTE THE KOCHIN FUNCTIONS
theta_range = np.linspace(0, 2*pi, 40)
omega_range = np.linspace(1.0, 5.0, 3)

test_matrix = xr.Dataset(coords={
    'omega': omega_range, 'theta': theta_range, 'radiating_dof': ["Surge"],
})
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, buoy.immersed_part(), wavenumber=True)

# PLOT THE KOCHIN FUNCTION
plt.figure()
dataset = dataset.swap_dims({'omega': 'wavenumber'})  # Use wavenumber to index data
for k in dataset.coords['wavenumber']:
    (dataset['kochin']
     .sel(radiating_dof="Surge", wavenumber=k)
     .real
     .plot(x='theta', label=f"k={float(k):.2f} m¯¹", ax=plt.gca()))
plt.legend()
plt.title("Real part of the Kochin function as a function of theta")
plt.tight_layout()
plt.show()
