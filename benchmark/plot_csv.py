#!/usr/bin/env python
# coding: utf-8

import sys
import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

directory = sys.argv[1]
# directory = '2017-10-02_173227/'

#######################################################################
#                              Plot time                              #
#######################################################################

times = pd.DataFrame.from_csv(os.path.join(directory, 'times.csv'))

times = times.groupby('nb_cells').aggregate(np.mean)
ax = times.plot()
ax.set(xlabel='number of cells in mesh', ylabel='time (seconds)')

plt.grid()

#######################################################################
#                            Check values                             #
#######################################################################
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

# added_mass.plot(y=[name for name in case_names if '600' in name])
# print(added_mass[[name for name in case_names if '600' in name]])

plt.show()
