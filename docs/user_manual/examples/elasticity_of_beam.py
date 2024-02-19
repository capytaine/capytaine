#!/usr/bin/env python

import numpy as np
import xarray as xr
import capytaine as cpt
import matplotlib.pyplot as plt
from scipy.optimize import newton

# We consider a vertical cylinder clamped at the sea bottom and freely moving at the top above the free surface.

water_depth = 10

cylinder_length = 12.0
cylinder_density = 1000.0
cylinder_diameter = 1.0
young_modulus = 1e9
inertia_moment = np.pi * cylinder_diameter**4 / 64  # Cylinder with circular section

# The following code computes the shapes of the deformation modes of the body and the associated eigenfrequencies.

def deformation_mode(i_mode, z):
    alpha = newton(lambda a: np.cos(a)*np.cosh(a) + 1, (2*i_mode + 1)*np.pi/2)
    coeff = (-np.cos(alpha) - np.cosh(alpha))/(np.sin(alpha) + np.sinh(alpha))
    k = alpha/cylinder_length
    z_ = z + water_depth
    return 0.5 * (np.cos(k*z_) - np.cosh(k*z_) + coeff*(np.sin(k*z_) - np.sinh(k*z_)))

def eigenfrequency(i_mode):
    alpha = newton(lambda a: np.cos(a)*np.cosh(a) + 1, (2*i_mode + 1)*np.pi/2)
    return alpha/cylinder_length**2 * np.sqrt(young_modulus*inertia_moment/(cylinder_density*np.pi*(cylinder_diameter/2)**2))


# Plotting the x-deformation of the first 4 modes
z_range = np.linspace(-water_depth, -water_depth + cylinder_length, 100)

# for i_mode in range(4):
#     plt.plot(deformation_mode(i_mode, z_range), z_range, label=f"mode_{i_mode}")
# plt.legend()


mesh = cpt.mesh_vertical_cylinder(
        center=(0.0, 0.0, cylinder_length/2 - water_depth - 0.1),
        length=cylinder_length,
        radius=cylinder_diameter/2,
        resolution=(4, 15, 36),
        )

nb_modes = 10
dofs = {}
for i in range(nb_modes):
    dofs[f"mode_{i}"] = np.array([(deformation_mode(i, z), 0.0, 0.0) for x, y, z in mesh.faces_centers])

full_body = cpt.FloatingBody(mesh=mesh, dofs=dofs, center_of_mass=(0, 0, cylinder_length/2 - water_depth))
body = full_body.immersed_part(water_depth=water_depth)

body.inertia_matrix = body.add_dofs_labels_to_matrix(np.eye(nb_modes))
body.internal_stiffness_matrix = body.add_dofs_labels_to_matrix(np.diag([eigenfrequency(i_mode) for i_mode in range(nb_modes)]))
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

test_matrix = xr.Dataset(coords={
    "omega": [1.0],
    "radiating_dof": list(body.dofs),
    "wave_direction": [0.0],
    "water_depth": [water_depth],
    })
ds = cpt.BEMSolver().fill_dataset(test_matrix, body)
rao = cpt.post_pro.rao(ds, stiffness=body.internal_stiffness_matrix)


deformation_profile = sum(
    rao.isel(radiating_dof=i_mode).values[0] * deformation_mode(i_mode, z_range)
    for i_mode, dof in enumerate(body.dofs)
    )
plt.plot(np.real(deformation_profile), z_range)

plt.show()

