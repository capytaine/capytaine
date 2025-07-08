import numpy as np
import xarray as xr
import capytaine as cpt
import matplotlib.pyplot as plt

cpt.set_logging('INFO')

# Initialize floating body by generating a geometric mesh
mesh = cpt.mesh_horizontal_cylinder(
    length=10.0, radius=0.5,  # Dimensions
    center=(0, 0, 0),         # Position
    resolution=(5, 20, 40)    # Fineness of the mesh
    ).immersed_part()
body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())

# Define the range of water depth
depth_range = np.concatenate([np.linspace(2, 100, 10), [np.inf]])

test_matrix = xr.Dataset(coords={
    "wavelength": 100.0,
    "water_depth": depth_range,
    "radiating_dof": ["Heave"]
    })

# Solve all radiation problems
solver = cpt.BEMSolver()
data = solver.fill_dataset(test_matrix, body)

# Note that the solver could not solve the most shallow case and returned NaN:
data.added_mass.sel(radiating_dof="Heave", influenced_dof="Heave", water_depth=2.0).values

# Plot the added mass of each dofs as a function of the water depth.
fig, ax = plt.subplots(layout="constrained")
ax.plot(
    data.water_depth * data.wavenumber,
    data['added_mass'].sel(radiating_dof="Heave", influenced_dof="Heave"),
    marker="s", label="Finite depth"
)
ax.hlines(
    data['added_mass'].sel(radiating_dof="Heave", influenced_dof="Heave", water_depth=np.inf),
    xmin=5, xmax=8, label="Infinite depth",
    linestyle="--"
)
ax.set(xlabel='wavenumber * water_depth', ylabel='Heave-heave added mass')
ax.legend()
plt.show()
