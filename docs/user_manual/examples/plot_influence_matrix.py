#!/usr/bin/env python

import numpy as np
import capytaine as cpt

# Generate the mesh of a cylinder
cylinder = cpt.HorizontalCylinder(
    length=10.0, radius=1.0,  # Dimensions
    center=(0, 0, -2),        # Position
    nr=1, nx=8, ntheta=6,     # Fineness of the mesh
)

engine = cpt.BasicMatrixEngine()
green_function = cpt.Delhommeau()

S, K = engine.build_matrices(
    cylinder.mesh, cylinder.mesh,
    free_surface=0.0, sea_bottom=-np.infty,
    wavenumber=1.0,
    green_function=green_function,
)

# Plot the absolute value of the matrix S
import matplotlib.pyplot as plt
plt.imshow(abs(S.full_matrix()))
plt.colorbar()
plt.title("$|S|$")
plt.tight_layout()
plt.show()
