#!/usr/bin/env python

import logging

import numpy as np

from capytaine import (FloatingBody, AxialSymmetry, Nemoh,
                       RadiationProblem, assemble_dataset)

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

# Profile of the axisymmetric body
def shape(z):
    return 0.1*(-(z+1)**2 + 16)

# Generate the mesh and display it with VTK.
buoy = FloatingBody(
    AxialSymmetry.from_profile(shape, z_range=np.linspace(-5, 0, 30), nphi=40)
)
buoy.add_translation_dof(name="Heave")
buoy.show()

# Set up problems
omega_range = np.linspace(0.1, 5.0, 60)
problems = [RadiationProblem(body=buoy, radiating_dof='Heave', omega=omega)
            for omega in omega_range]

# Solve the problems
solver = Nemoh()
results = [solver.solve(pb) for pb in sorted(problems)]
dataset = assemble_dataset(results)

# Plot results
import matplotlib.pyplot as plt
plt.figure()
plt.plot(
    omega_range,
    dataset['added_mass'].sel(radiating_dof='Heave',
                              influenced_dof='Heave'),
    label="Added mass",
)
plt.plot(
    omega_range,
    dataset['radiation_damping'].sel(radiating_dof='Heave',
                                     influenced_dof='Heave'),
    label="Radiation damping",
)
plt.xlabel('omega')
plt.grid()
plt.legend()
plt.show()
