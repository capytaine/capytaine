#!/usr/bin/env python
# coding: utf-8
"""
"""

import numpy as np
import matplotlib.pyplot as plt

from capytaine.reference_bodies import *
from capytaine.symmetries import *
from capytaine.Nemoh import Nemoh

# cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=10, nr=0, ntheta=10)
# cylinder.translate_z(-2.0)
# cylinder.show()

nb_slices = 5
nb_theta = 10
ring = HorizontalCylinder(length=10.0/nb_slices, radius=1.0, nx=2, nr=0, ntheta=nb_theta+1)
ring.translate_x(-5.0)
ring.translate_z(-2.0)
ring.name = "ring"
cylinder = TranslationalSymmetry(ring, translation=np.asarray([10.0/nb_slices, 0.0, 0.0]), nb_repetitions=nb_slices-1)
cylinder.show()

solver = Nemoh()

S, V = cylinder.build_matrices(cylinder)

plt.imshow(np.abs(V.full_matrix()))
plt.colorbar()
plt.title("$|V|$")
plt.tight_layout()
plt.show()
