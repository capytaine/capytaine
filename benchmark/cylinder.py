#!/usr/bin/env python
# coding: utf-8

import datetime

import numpy as np
import pandas as pd

from capytaine.reference_bodies import HorizontalCylinder
from capytaine.symmetries import *

from benchmark import *

WORKING_DIRECTORY = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

def full_resolution_Nemoh(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=nb_slices+1, nr=0, ntheta=nb_theta+1)
    cylinder.translate_x(-5.0)
    cylinder.translate_z(-2.0)

    return profile_Nemoh(cylinder, omega_range, f"{WORKING_DIRECTORY}/Nemoh_{nb_theta*nb_slices}")

def full_resolution_Capytaine(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, nx=nb_slices+1, nr=0, ntheta=nb_theta+1)
    cylinder.translate_x(-5.0)
    cylinder.translate_z(-2.0)
    cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)

    return profile_capytaine(cylinder, omega_range, f"{WORKING_DIRECTORY}/Capy_{nb_theta*nb_slices}")

def clever_resolution_Capytaine(nb_slices, nb_theta, omega_range):
    ring = HorizontalCylinder(length=10.0/nb_slices, radius=1.0, nx=2, nr=0, ntheta=nb_theta+1)
    ring.translate_x(-5.0)
    ring.translate_z(-2.0)
    trans_cylinder = TranslationSymmetry(ring, translation=np.asarray([10.0/nb_slices, 0.0, 0.0]), nb_repetitions=nb_slices-1)
    trans_cylinder.dofs["Heave"] = trans_cylinder.faces_normals @ (0, 0, 1)

    return profile_capytaine(trans_cylinder, omega_range, f"{WORKING_DIRECTORY}/Capy_{nb_theta*nb_slices}")


nb_theta = 10
nb_slices_range = range(10, 100, 5)
omega_range = np.linspace(0.1, 4.0, 40)
times = pd.DataFrame(
    index=(nb_slices*nb_theta for nb_slices in nb_slices_range),
    columns=["Nemoh 2.0", "Refactored", "Refactored + symmetries"]
)

for nb_slices in nb_slices_range:
    print(f"====> Slices: {nb_slices}")
    times["Nemoh 2.0"][nb_slices*nb_theta] = full_resolution_Nemoh(nb_slices, nb_theta, omega_range)
    times["Refactored"][nb_slices*nb_theta] = full_resolution_Capytaine(nb_slices, nb_theta, omega_range)
    times["Refactored + symmetries"][nb_slices*nb_theta] = clever_resolution_Capytaine(nb_slices, nb_theta, omega_range)

print(times)
times.to_csv(f'{WORKING_DIRECTORY}/times.csv')
