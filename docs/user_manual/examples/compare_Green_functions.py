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
    'omega': np.linspace(0.5, 4, 100),
    'radiating_dof': list(body.dofs.keys()),
})

green_functions = [
        cpt.Delhommeau(gf_singularities="low_freq"),
        # cpt.Delhommeau(gf_singularities="high_freq"),  # For this problem, more difficult to converge
        cpt.LiangWuNoblesseGF(),
        ]

data = []
for gf in green_functions:
    data.append(cpt.BEMSolver(green_function=gf).fill_dataset(test_matrix, body))

fig, axs = plt.subplots(2, 1, sharex=True, layout="constrained")
for gf, ds in zip(green_functions, data):
    ds['added_mass'].plot(ax=axs[0], x='omega', label=str(gf))
    ds['radiation_damping'].plot(ax=axs[1], x='omega', label=str(gf))
axs[0].set_title("Added mass")
axs[0].legend()
axs[1].set_title("Radiation damping")

plt.show()
