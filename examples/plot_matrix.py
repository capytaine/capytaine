#!/usr/bin/env python
# coding: utf-8
"""
"""

import numpy as np
import matplotlib.pyplot as plt

from capytaine import *

cylinder = generate_clever_horizontal_cylinder(length=10.0, radius=1.0, nx=5, ntheta=10)
cylinder.translate_z(-2.0)
# cylinder.show()

solver = Nemoh()

S, V = solver.build_matrices(cylinder, cylinder, free_surface=np.infty)
# S, V = cylinder.build_matrices(cylinder)

plt.imshow(np.abs(V.full_matrix()))
plt.colorbar()
plt.title("$|V|$")
plt.tight_layout()
plt.show()
