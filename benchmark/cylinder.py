#!/usr/bin/env python
# coding: utf-8

import datetime
from itertools import count

import numpy as np
import pandas as pd

from capytaine.reference_bodies import *
from capytaine.symmetries import *

from benchmark import profile_Nemoh, profile_capytaine

WORKING_DIRECTORY = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

ID = count()

PROBLEM_ARGS = {
    "sea_bottom": -6.0,
    "rho": 1000,
}

CYLINDER_POSITION = np.array((-5.0, 0.0, -3.0))

def full_resolution_Nemoh(nb_slices, nb_theta, omega_range):
    cylinder = generate_open_horizontal_cylinder(length=10.0, radius=1.0, nx=nb_slices, ntheta=nb_theta)
    cylinder.translate(CYLINDER_POSITION)
    return profile_Nemoh(cylinder,
                         omega_range,
                         f"{WORKING_DIRECTORY}/{next(ID):03}_Nemoh_{nb_theta*nb_slices}",
                         **PROBLEM_ARGS)


def full_Capytaine(nb_slices, nb_theta, omega_range):
    cylinder = generate_open_horizontal_cylinder(length=10.0, radius=1.0, nx=nb_slices, ntheta=nb_theta)
    cylinder.translate(CYLINDER_POSITION)
    cylinder.dofs["Heave"] = cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_Capytaine(nb_slices, nb_theta, omega_range):
    half_cylinder = generate_open_horizontal_cylinder(length=10.0, radius=1.0, nx=nb_slices, ntheta=nb_theta//2, half=True)
    half_cylinder.translate(CYLINDER_POSITION)
    sym_cylinder = ReflectionSymmetry(half_cylinder, xOz_Plane)
    sym_cylinder.dofs["Heave"] = sym_cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(sym_cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_sym_Capytaine(nb_slices, nb_theta, omega_range):
    half_cylinder = generate_open_horizontal_cylinder(length=10.0, radius=1.0, nx=nb_slices, ntheta=nb_theta//2, half=True)
    half_cylinder.translate(CYLINDER_POSITION)
    quarter_cylinder = half_cylinder.extract_faces(np.where(half_cylinder.faces_centers[:, 0] > 0)[0]) # Keep x > 0
    quarter_cylinder.name = "quarter_cylinder"
    sym_sym_cylinder = ReflectionSymmetry(ReflectionSymmetry(quarter_cylinder, xOz_Plane), yOz_Plane)
    sym_sym_cylinder.dofs["Heave"] = sym_sym_cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(sym_sym_cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_sym_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def trans_Capytaine(nb_slices, nb_theta, omega_range):
    ring = generate_ring(length=10.0/nb_slices, radius=1.0, ntheta=nb_theta)
    trans_cylinder = TranslationalSymmetry(ring, translation=np.asarray([10.0/nb_slices, 0.0, 0.0]), nb_repetitions=nb_slices-1)
    trans_cylinder.translate(CYLINDER_POSITION)
    trans_cylinder.dofs["Heave"] = trans_cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(trans_cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_trans_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def trans_sym_Capytaine(nb_slices, nb_theta, omega_range):
    half_ring = generate_ring(length=10.0/nb_slices, radius=1.0, ntheta=nb_theta//2, half=True)
    half_ring.translate(CYLINDER_POSITION)
    trans_cylinder = ReflectionSymmetry(TranslationalSymmetry(half_ring, translation=np.asarray([10.0/nb_slices, 0.0, 0.0]), nb_repetitions=nb_slices-1), xOz_Plane)
    trans_cylinder.dofs["Heave"] = trans_cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(trans_cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_trans_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_trans_Capytaine(nb_slices, nb_theta, omega_range):
    half_ring = generate_ring(length=10.0/nb_slices, radius=1.0, ntheta=nb_theta//2, half=True)
    half_ring.translate(CYLINDER_POSITION)
    trans_cylinder = TranslationalSymmetry(ReflectionSymmetry(half_ring, xOz_Plane), translation=np.asarray([10.0/nb_slices, 0.0, 0.0]), nb_repetitions=nb_slices-1)
    trans_cylinder.dofs["Heave"] = trans_cylinder.faces_normals @ (0, 0, 1)
    return profile_capytaine(trans_cylinder,
                             omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_sym_trans_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


# ===============================================================
# ===============================================================
# ===============================================================

if __name__ == "__main__":

    nb_repetitions = 3
    nb_cells_range = 2*(np.sqrt(np.linspace(200, 3000, 15))//2)
    nb_theta_range = [int(x) for x in nb_cells_range]
    nb_slices_range = [int(x) for x in nb_cells_range]
    omega_range = np.linspace(0.1, 4.0, 40)


    cases = {
        # "Nemoh 2.0":                          full_resolution_Nemoh,
        "Capytaine":                          full_Capytaine,
        # "Capytaine + symmetry":               sym_Capytaine,
        # "Capytaine + 2 symmetries":           sym_sym_Capytaine,
        "Capytaine + translation":            trans_Capytaine,
        # "Capytaine + translation + symmetry": trans_sym_Capytaine,
        # "Capytaine + symmetry + translation": sym_trans_Capytaine,
    }

    # ===========================

    times = pd.DataFrame(
        index=range(nb_repetitions*len(nb_slices_range)),
        columns=["nb_cells"] + list(cases.keys()),
    )

    for i, ((nb_slices, nb_theta), k) in enumerate(product(zip(nb_slices_range, nb_theta_range), range(nb_repetitions))):
        print(f"====> {nb_slices}×{nb_theta} cells, Repetition: {k+1}/{nb_repetitions}")
        times["nb_cells"][i] = nb_slices*nb_theta
        for name, function in cases.items():
            print("\t\t" + name)
            times[name][i] = function(nb_slices, nb_theta, omega_range)

    print(times)
    times.to_csv(f'{WORKING_DIRECTORY}/times.csv')
