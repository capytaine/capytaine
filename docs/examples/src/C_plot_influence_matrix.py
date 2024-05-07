#!/usr/bin/env python

import numpy as np
import capytaine as cpt

# Generate the mesh of a cylinder
cylinder = cpt.FloatingBody(
        mesh=cpt.mesh_horizontal_cylinder(
            length=10.0, radius=1.0,
            center=(0, 0, -2),
            resolution=(1, 6, 8)
            ))

engine = cpt.BasicMatrixEngine()
green_function = cpt.Delhommeau()

S, K = engine.build_matrices(
    cylinder.mesh, cylinder.mesh,
    free_surface=0.0, water_depth=np.inf,
    wavenumber=1.0,
    green_function=green_function,
)

# Plot the absolute value of the matrix S
import matplotlib.pyplot as plt
plt.imshow(abs(np.array(S)))
plt.colorbar()
plt.title("$|S|$")
plt.tight_layout()
plt.show()
