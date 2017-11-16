#!/usr/bin/env python
# coding: utf-8

import sys
import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#######################################################################
#                              Plot time                              #
#######################################################################
def compare_all_total_times(directory):
    times = pd.DataFrame.from_csv(os.path.join(directory, 'times.csv'))

    # print(times.groupby('nb_cells').aggregate(np.std))

    times = times.groupby('nb_cells').aggregate(np.min)
    ax = times.plot()
    ax.set(
        xlabel='number of cells in mesh',
        ylabel='computation time (seconds)',
    )
    plt.grid()
    plt.tight_layout()


#######################################################################
#                            Time details                             #
#######################################################################
def plot_detailed_time(directory):
    dirs = glob.glob(os.path.join(directory, '*capy*'))
    detailed_time = pd.DataFrame(
        index=dirs,
        columns=['type', 'nb_cells', 'evaluate matrices', 'solve linear problem', 'total'],
    )

    for result_dir in dirs:

        detailed_time['nb_cells'][result_dir] = int(result_dir.split('_')[-1])
        detailed_time['type'][result_dir] = '_'.join(result_dir.split('_')[2:-1])

        with open(os.path.join(result_dir, 'profile.log'), 'r') as profile_file:
            for entry in profile_file.readlines():
                if '(build_matrices)' in entry:
                    detailed_time['evaluate matrices'][result_dir] = float(entry.split()[3])
                elif '(solve)' in entry and 'numpy/linalg' in entry:
                    detailed_time['solve linear problem'][result_dir] = float(entry.split()[3])
                elif '(solve)' in entry and 'Nemoh.py' in entry:
                    detailed_time['total'][result_dir] = float(entry.split()[3])

    detailed_time['other'] = detailed_time['total'] - detailed_time['evaluate matrices'] - detailed_time['solve linear problem']
    detailed_time = detailed_time.sort_values(by='nb_cells')
    detailed_time = detailed_time.groupby(['type', 'nb_cells']).aggregate(np.min)

    for test_type in detailed_time.index.levels[0]:
        ax = detailed_time.T[test_type].T.plot.area(y=['solve linear problem', 'evaluate matrices', 'other'])
        ax.set(
            ylim=(0.0, 80.0),
            xlabel='number of cells in mesh',
            ylabel='computation time (seconds)',
        )
        plt.tight_layout()
        plt.grid(zorder=3)


#######################################################################
#                            Check values                             #
#######################################################################
def compare_results(directory):
    omega_range = np.linspace(0.1, 4.0, 40)

    nemoh_dirs = glob.glob(os.path.join(directory, '*Nemoh*'))
    capy_dirs = glob.glob(os.path.join(directory, '*capy*'))

    case_names = sorted(os.path.basename(name) for name in nemoh_dirs+capy_dirs)

    added_mass = pd.DataFrame(index=omega_range, columns=case_names)
    damping = pd.DataFrame(index=omega_range, columns=case_names)

    for nemoh_dir in nemoh_dirs:
        mesh = int(nemoh_dir.split('_')[-1])
        results = np.genfromtxt(os.path.join(nemoh_dir, "results", "Forces.dat"))
        added_mass[os.path.basename(nemoh_dir)] = results[::2]
        damping[os.path.basename(nemoh_dir)] = results[1::2]

    for capy_dir in capy_dirs:
        mesh = int(capy_dir.split('_')[-1])
        results = np.genfromtxt(os.path.join(capy_dir, "results.csv"))
        added_mass[os.path.basename(capy_dir)] = results[:, 0]
        damping[os.path.basename(capy_dir)] = results[:, 1]

    added_mass.plot(y=[name for name in case_names if '600' in name])
    # print(added_mass[[name for name in case_names if '600' in name]])


if __name__ == "__main__":
    plot_detailed_time(sys.argv[1])
    compare_all_total_times(sys.argv[1])
    compare_results(sys.argv[1])

    # plot_detailed_time("2017-11-06_153544/")
    # compare_all_total_times("2017-11-06_153544/")

    plt.show()

