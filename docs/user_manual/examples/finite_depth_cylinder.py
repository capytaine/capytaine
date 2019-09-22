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

# Define a degree of freedom. The keyword "Heave"
# is recognized by the code and the vertical translation
# motion is automatically defined.
cylinder.add_translation_dof(name="Heave")

# Define the range of water depth
depth_range = list(range(5, 25, 2)) + [np.infty]

# Set up the problems: we will solve a radiation problem for each
# water depth:
problems = [
    cpt.RadiationProblem(body=cylinder, sea_bottom=-depth, omega=2.0)
    for depth in depth_range
]
# Water density, gravity and radiating dof have not been specified.
# Default values are used. (For the radiating dof, the default value
# is usually the first one that has been defined. Here only one has
# been defined.)

# Solve all radiation problems
solver = cpt.BEMSolver(engine=cpt.HierarchicalToeplitzMatrixEngine())
results = [solver.solve(pb) for pb in sorted(problems)]

# Gather the computed added mass into a labelled array.
data = cpt.assemble_dataset(results)

# Plot the added mass of each dofs as a function of the water depth.
import matplotlib.pyplot as plt
plt.figure()
plt.plot(
    depth_range,
    data['added_mass'].sel(omega=2.0, radiating_dof="Heave", influenced_dof="Heave"),
    marker="s",
)
plt.xlabel('water depth')
plt.ylabel('added mass')
plt.show()
