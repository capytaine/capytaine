#!/usr/bin/env python

import numpy as np

from capytaine import Nemoh
from capytaine.geometric_bodies import HorizontalCylinder

# Generate the mesh of a cylinder
cylinder = HorizontalCylinder(
    length=10.0, radius=1.0,  # Dimensions
    center=(0, 0, -2),        # Position
    nr=1, nx=8, ntheta=6,     # Fineness of the mesh
)

# Use Nemoh to compute the influence matrices
solver = Nemoh()
S, K = solver.build_matrices(cylinder.mesh, cylinder.mesh, wavenumber=1.0)

# Plot the absolute value of the matrix V
#
import matplotlib.pyplot as plt
plt.imshow(np.abs(S.full_matrix()))
plt.colorbar()
plt.title("$|S|$")
plt.tight_layout()
plt.show()
