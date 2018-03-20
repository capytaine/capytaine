#!/usr/bin/env python
# coding: utf-8
"""
Example computation: added mass and damping of an oscillating flap.
"""

import os
import glob

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from capytaine import *

result_directory = os.path.join(os.path.dirname(__file__), "flap_results")


def solve_flap(clever=True, resolution=2):
    """Solve the flap problem for a given resolution.

    Parameters
    ----------
    clever : bool, optional
        if True, use the prismatic symmetry of the flap to speed up the computations
    resolution: int, optional
        the number of cells in the mesh will be proportional to this coefficient
    """
    # Load reference range of frequencies
    T_range, _, _ = np.loadtxt(os.path.join(os.path.dirname(__file__), "data/flap_mu_nu.tsv")).T

    depth = 10.9

    # Create mesh
    if clever:
        # Use prismatic shape to speed up computations.
        flap = generate_clever_open_rectangular_parallelepiped(
            height=depth, width=3.0, thickness=0.001,
            nh=int(10*resolution), nw=int(3*resolution))
    else:
        # Do not use prismatic shape to speed up computations.
        flap = generate_open_rectangular_parallelepiped(
            height=depth, width=3.0, thickness=0.001,
            nh=int(10*resolution), nw=int(3*resolution))

    flap.translate_z(-depth)

    # Set oscillation degree of freedom
    # The upper part of the flap rotation around an horizontal axis.
    # The lower part is fixed.
    flap.add_rotation_dof(axis_direction=(0, 1, 0), axis_point=(0, 0, -9.4), name="Oscillation")
    flap.dofs["Oscillation"][flap.mesh.faces_centers[:, 2] < -9.4] = 0.0

    # Set up problems and initialise solver
    problems = [RadiationProblem(body=flap, omega=omega, radiating_dof="Oscillation", sea_bottom=-depth) for omega in 2*np.pi/T_range]
    solver = Nemoh()

    # Solve problems
    results = [solver.solve(problem) for problem in problems]
    return assemble_dataset(results)


def plot_flap_results():
    """Plot the results of all the csv files in result directory."""

    # Load reference result
    T_range, mu, nu = np.loadtxt(os.path.join(os.path.dirname(__file__), "data/flap_mu_nu.tsv")).T

    # Plot reference results
    plt.figure(1)
    plt.plot(T_range, mu, linestyle="--", label="Reference added mass")

    plt.figure(2)
    plt.plot(T_range, nu, linestyle="--", label="Reference added damping")

    # Plot all results in result directory
    for fichier in glob.glob(os.path.join(result_directory, "*.nc")):
        dataset = xr.open_dataset(fichier)

        nb_cells = int(fichier.split('_')[-2])

        plt.figure(1)
        plt.plot(T_range, dataset['added_mass'][::-1, 0, 0], label=f"Added mass ({nb_cells} cells)")

        plt.figure(2)
        plt.plot(T_range, dataset['radiation_damping'][::-1, 0, 0], label=f"Added damping ({nb_cells} cells)")

    plt.figure(1)
    plt.xlabel("Wave period (s)")
    plt.legend()
    plt.tight_layout()

    plt.figure(2)
    plt.xlabel("Wave period (s)")
    plt.legend()
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    resolution = 2
    clever = True

    import logging
    from datetime import datetime

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s:\t%(message)s", datefmt="%H:%M:%S")
    LOG = logging.getLogger(__name__)

    # Create directory to store results
    if not os.path.isdir(result_directory):
        os.mkdir(result_directory)
    result_file_path = os.path.join(os.path.dirname(__file__), "flap_results",
                                    f"Results_{'clever_' if clever else ''}{30*resolution**2}_cells.nc")
    if os.path.exists(result_file_path):
        LOG.warning(f"Overwriting {result_file_path}")

    start_time = datetime.now()
    dataset = solve_flap(resolution, clever)
    dataset.to_netcdf(result_file_path)
    end_time = datetime.now()
    LOG.info(f"Duration: {end_time - start_time}")

    plot_flap_results()

