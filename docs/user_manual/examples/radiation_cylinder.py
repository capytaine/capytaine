#!/usr/bin/env python

import logging
import numpy as np
import capytaine as cpt

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

# Initialize floating body by generating a geometric mesh
cylinder = cpt.HorizontalCylinder(
    length=10.0, radius=1.0,  # Dimensions
    center=(0, 0, -2),        # Position
    nr=5, nx=40, ntheta=20,   # Fineness of the mesh
)

# Automatically add the six degrees of freedom of a rigid body
cylinder.add_all_rigid_body_dofs()

# Define the range of frequencies as a Numpy array
omega_range = np.linspace(0.1, 6.0, 20)

# Set up the problems: we will solve a radiation problem for each
# degree of freedom of the body and for each frequency in the
# frequency range.
problems = [
    cpt.RadiationProblem(body=cylinder, radiating_dof=dof, omega=omega)
    for dof in cylinder.dofs
    for omega in omega_range
]
# Water density, gravity and water depth have not been specified.
# Default values are used.

# Solve all radiation problems
solver = cpt.BEMSolver()
results = [solver.solve(pb) for pb in sorted(problems)]
# The 'sorted' function ensures that the problems are sequentially
# treated in an optimal order.

# Gather the computed added mass into a labelled array.
data = cpt.assemble_dataset(results)

# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
for dof in cylinder.dofs:
    plt.plot(
        omega_range,
        data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )
plt.xlabel('omega')
plt.ylabel('added mass')
plt.legend()
plt.tight_layout()
plt.show()
