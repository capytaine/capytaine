#!/usr/bin/env python
# coding: utf-8

import datetime
from itertools import count, product

import numpy as np
import pandas as pd

from capytaine.reference_bodies import generate_axi_symmetric_body

from benchmark import *

WORKING_DIRECTORY = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

ID = count()


def shape(z):
    return 0.1*(-(z+1)**2 + 16)


def full_resolution_Nemoh(nb_slices, nz, omega_range):
    buoy = generate_axi_symmetric_body(shape, z_range=np.linspace(-5.0, 0.0, nz), nphi=nb_slices)
    buoy.dofs["Heave"] = buoy.faces_normals @ (0, 0, 1)
    return profile_Nemoh(buoy.as_FloatingBody(), omega_range,
                         f"{WORKING_DIRECTORY}/{next(ID):03}_Nemoh_{nz*nb_slices}")


def full_Capytaine(nb_slices, nz, omega_range):
    buoy = generate_axi_symmetric_body(shape, z_range=np.linspace(-5.0, 0.0, nz), nphi=nb_slices)
    buoy.dofs["Heave"] = buoy.faces_normals @ (0, 0, 1)
    return profile_capytaine(buoy.as_FloatingBody(), omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_capy_{nz*nb_slices}")


def rot_Capytaine(nb_slices, nz, omega_range):
    buoy = generate_axi_symmetric_body(shape, z_range=np.linspace(-5.0, 0.0, nz), nphi=nb_slices)
    buoy.dofs["Heave"] = buoy.faces_normals @ (0, 0, 1)
    return profile_capytaine(buoy, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_rot_capy_{nz*nb_slices}")


# ===============================================================
# ===============================================================
# ===============================================================

if __name__ == "__main__":

    nb_repetitions = 1
    nz = 30
    nb_slices_range = range(60, 90, 10)
    omega_range = np.linspace(0.1, 4.0, 40)

    cases = {
        "Nemoh 2.0":                full_resolution_Nemoh,
        "Capytaine":                full_Capytaine,
        "Capytaine + axisymmetry":  rot_Capytaine,
        }

    # ===========================

    times = pd.DataFrame(
        index=range(nb_repetitions*len(nb_slices_range)),
        columns=["nb_cells"] + list(cases.keys()),
    )

    for i, (nb_slices, k) in enumerate(product(nb_slices_range, range(nb_repetitions))):
        print(f"====> Slices: {nb_slices}, Repetition: {k+1}/{nb_repetitions}")
        times["nb_cells"][i] = nb_slices*nz
        for name, function in cases.items():
            print("\t\t" + name)
            times[name][i] = function(nb_slices, nz, omega_range)

    print(times)
    times.to_csv(f'{WORKING_DIRECTORY}/times.csv')
