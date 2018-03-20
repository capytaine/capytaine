#!/usr/bin/env python
# coding: utf-8

import datetime
from itertools import count

import pandas as pd

from capytaine.bodies import FloatingBody
from capytaine.geometric_bodies.cylinder import HorizontalCylinder
from capytaine.symmetries import ReflectionSymmetry

from benchmark import profile_Nemoh, profile_capytaine

WORKING_DIRECTORY = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

ID = count()

PROBLEM_ARGS = {
    "sea_bottom": -6.0,
    "rho": 1000,
}

CYLINDER_POSITION = np.array((-0.0, 0.0, -3.0))

def full_resolution_Nemoh(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, center=CYLINDER_POSITION,
                                  clever=False, nx=nb_slices, nr=0, ntheta=nb_theta)
    cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_Nemoh(cylinder, omega_range,
                         f"{WORKING_DIRECTORY}/{next(ID):03}_Nemoh_{nb_theta*nb_slices}",
                         **PROBLEM_ARGS)


def full_Capytaine(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, center=CYLINDER_POSITION,
                                  clever=False, nx=nb_slices, nr=0, ntheta=nb_theta)
    cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(cylinder, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_Capytaine(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, center=CYLINDER_POSITION,
                                  clever=False, nx=nb_slices, nr=0, ntheta=nb_theta)
    half_cylinder_mesh = cylinder.mesh.extract_faces(np.where(cylinder.mesh.faces_centers[:, 1] > 0)[0]) # Keep y > 0
    half_cylinder.name = "half_cylinder"
    sym_cylinder = FloatingBody(ReflectionSymmetry(half_cylinder_mesh, plane=xOz_Plane))
    sym_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(sym_cylinder, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_sym_Capytaine(nb_slices, nb_theta, omega_range):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, center=CYLINDER_POSITION,
                                  clever=False, nx=nb_slices, nr=0, ntheta=nb_theta)
    half_cylinder_mesh = cylinder.mesh.extract_faces(np.where(cylinder.mesh.faces_centers[:, 1] > 0)[0]) # Keep y > 0
    quarter_cylinder_mesh = half_cylinder_mesh.extract_faces(np.where(half_cylinder_mesh.faces_centers[:, 0] > 0)[0]) # Keep x > 0
    quarter_cylinder.name = "quarter_cylinder"
    sym_sym_cylinder = ReflectionSymmetry(ReflectionSymmetry(quarter_cylinder, plane=xOz_Plane), plane=yOz_Plane)
    sym_sym_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(sym_sym_cylinder, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_sym_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def trans_Capytaine(nb_slices, nb_theta, omega_range):
    ring = HorizontalCylinder(length=10.0/nb_slices, radius=1.0, center=CYLINDER_POSITION - np.array((5.0, 0, 0),
                                  clever=False, nx=1, nr=0, ntheta=nb_theta)
    trans_cylinder = FloatingBody(TranslationalSymmetry(ring, translation=(10.0/nb_slices, 0.0, 0.0), nb_repetitions=nb_slices-1))
    trans_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(trans_cylinder, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_trans_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def trans_sym_Capytaine(nb_slices, nb_theta, omega_range):
    ring = HorizontalCylinder(length=10.0/nb_slices, radius=1.0, center=CYLINDER_POSITION - np.array((5.0, 0, 0),
                                  clever=False, nx=1, nr=0, ntheta=nb_theta)
    half_ring_mesh = ring.mesh.extract_faces(np.where(ring.mesh.faces_centers[:, 1] > 0)[0]) # Keep y > 0
    half_cylinder_mesh = TranslationalSymmetry(half_ring_mesh, translation=(10.0/nb_slices, 0.0, 0.0), nb_repetitions=nb_slices-1)
    trans_cylinder_mesh = ReflectionSymmetry(half_cylinder_mesh, xOz_Plane)
    trans_cylinder = FloatingBody(trans_cylinder_mesh)
    trans_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(trans_cylinder, omega_range,
                             f"{WORKING_DIRECTORY}/{next(ID):03}_trans_sym_capy_{nb_theta*nb_slices}",
                             **PROBLEM_ARGS)


def sym_trans_Capytaine(nb_slices, nb_theta, omega_range):
    ring = HorizontalCylinder(length=10.0/nb_slices, radius=1.0, center=CYLINDER_POSITION - np.array((5.0, 0, 0),
                                  clever=False, nx=1, nr=0, ntheta=nb_theta)
    half_ring_mesh = ring.mesh.extract_faces(np.where(ring.mesh.faces_centers[:, 1] > 0)[0]) # Keep y > 0
    ring_mesh = ReflectionSymmetry(half_ring_mesh, xOz_Plane)
    trans_cylinder_mesh = TranslationalSymmetry(ring_mesh, translation=(10.0/nb_slices, 0.0, 0.0), nb_repetitions=nb_slices-1)
    trans_cylinder = FloatingBody(trans_cylinder_mesh)
    trans_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return profile_capytaine(trans_cylinder, omega_range,
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
