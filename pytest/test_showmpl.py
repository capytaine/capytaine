#!/usr/bin/env python
# coding: utf-8
"""Test of method show_matplotlib."""

import logging
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

# Initialize floating body by generating a geometric mesh
cylinder = cpt.HorizontalCylinder(
    length=10.0, radius=1.0,  # Dimensions
    center=(0, 0, -2),        # Position
    nr=5, nx=40, ntheta=20,   # Fineness of the mesh
)


# Define the frequency
omega = 1.

# Set up the problem
problem = cpt.DiffractionProblem(body=cylinder, wave_direction=0, omega=omega)
# Water density, gravity and water depth have not been specified.
# Default values are used.

# Solve the problem
solver = cpt.BEMSolver()
results = solver.solve(problem)

# Plot just the mesh
cylinder.show_matplotlib()

# Plot the mesh with the potential field, by specifying ax, a label
# and a colormap
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

import matplotlib.cm as cm
colormap = cm.get_cmap('cividis')

cylinder.show_matplotlib(color_field=np.abs(results.potential),
                        ax=ax, cbar_label='$\phi (m^2/s)$')

ax.set_title('Potential distribution')
plt.show()

# As above, but avoid specifying ax, use default settings
cylinder.show_matplotlib(color_field=np.abs(results.potential))
