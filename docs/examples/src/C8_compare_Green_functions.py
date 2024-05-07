import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

cpt.set_logging('INFO')

# Generate body
mesh = cpt.mesh_horizontal_cylinder(
    length=3.0, radius=1.0, center=(0, 0, -1.01), resolution=(5, 30, 15)
)
body = cpt.FloatingBody(mesh)
body.add_translation_dof(name="Heave")

test_matrix = xr.Dataset(coords={
    'omega': np.linspace(0.5, 4, 40),
    'radiating_dof': list(body.dofs.keys()),
})

ds2 = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="high_freq")).fill_dataset(test_matrix, body)
ds1 = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).fill_dataset(test_matrix, body)

plt.figure()
ds1['added_mass'].plot(x='omega', label='High freq singularities')
ds2['added_mass'].plot(x='omega', label='Low freq singularities')
plt.legend()

plt.figure()
ds1['radiation_damping'].plot(x='omega', label='High freq singularities')
ds2['radiation_damping'].plot(x='omega', label='Low freq singularities')
plt.legend()

plt.show()
