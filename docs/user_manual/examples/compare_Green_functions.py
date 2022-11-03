import logging
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

# Generate body
body = cpt.HorizontalCylinder(
    length=3.0, radius=1.0,  # Dimensions
    center=(0, 0, -1.01),     # Position
    nr=5, nx=15, ntheta=30,   # Fineness of the mesh
)
body.add_translation_dof(name="Heave")

test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.5, 4, 40),
    'radiating_dof': list(body.dofs.keys()),
})

ds2 = cpt.BEMSolver(green_function=cpt.XieDelhommeau()).fill_dataset(test_matrix, body)
ds1 = cpt.BEMSolver(green_function=cpt.Delhommeau()).fill_dataset(test_matrix, body)

plt.figure()
ds1['added_mass'].plot(x='omega', label='Delhommeau')
ds2['added_mass'].plot(x='omega', label='XieDelhommeau')
plt.legend()

plt.figure()
ds1['radiation_damping'].plot(x='omega', label='Delhommeau')
ds2['radiation_damping'].plot(x='omega', label='XieDelhommeau')
plt.legend()

plt.show()
