"""Single body with a rotation center different from the center of mass.
Also a second body that is just a fixed obstacle.
"""

import numpy as np
import capytaine as cpt
import xarray as xr
import matplotlib.pyplot as plt

cpt.set_logging('INFO')

water_depth = 10.0
flap_height = 9.0
flap_draft = 8.0
flap_width = 6.0
base_height = water_depth - flap_draft - 0.1

flap_mesh = cpt.mesh_parallelepiped(
        size=(flap_width, 0.5, flap_height),
        center=(0.0, 0.0, -flap_draft + flap_height/2),
        faces_max_radius=0.2,
        missing_sides={"top"}
        )
flap_body = cpt.FloatingBody(
    mesh=flap_mesh,
    lid_mesh=flap_mesh.generate_lid(),
    dofs=cpt.rigid_body_dofs(only=['Pitch'], rotation_center=(0, 0, -flap_draft)),
    # The Pitch dof is defined around the hinge at z=-flap_draft.
    center_of_mass=(0, 0, -flap_draft + flap_height/2),
    name="flap"
)

base_mesh = cpt.mesh_parallelepiped(
        size=(1.1*flap_width, 0.8, base_height),
        center=(0.0, 0.0, -water_depth + base_height/2),
        faces_max_radius=0.2,
        missing_sides={"bottom"}
        )
base_body = cpt.FloatingBody(
    mesh=base_mesh,
    # no lid for fully immersed base
    # no dofs for fixed body
    center_of_mass=(0, 0, -10.90),
    name="base"
)

full_flap = flap_body + base_body

test_matrix = xr.Dataset(coords={
    "wavelength": np.linspace(5.0, 30.0, 10),
    "radiating_dof": ["flap__Pitch"],
    "wave_direction": [0.0],
    "water_depth": [water_depth],
    "rho": [1000.0],
    })

solver = cpt.BEMSolver()
dataset = solver.fill_dataset(
    test_matrix,
    full_flap.immersed_part(water_depth=water_depth)
)

dataset.added_mass.sel(
    radiating_dof="flap__Pitch",
    influenced_dof="flap__Pitch"
).plot(
    x="wavelength"
)
plt.show()
