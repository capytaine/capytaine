#!/usr/bin/env python
# coding: utf-8

import sys
import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

times = pd.DataFrame.from_csv(os.path.join(sys.argv[1], 'times.csv'))
ax = times.plot()
ax.set(xlabel='number of cells in mesh', ylabel='time (seconds)')
plt.grid()

# omega_range = np.linspace(0.1, 4.0, 40)

# mesh_to_plot = '600'
# plt.figure()
# for nemoh_dir in glob.glob(os.path.join(sys.argv[1], f'*Nemoh*{mesh_to_plot}*')):
#     mesh = int(nemoh_dir.split('_')[-1])
#     results = np.genfromtxt(os.path.join(nemoh_dir, "results", "Forces.dat"))
#     added_mass = results[::2]
#     damping = results[1::2]

#     plt.plot(omega_range, added_mass, label=nemoh_dir)

# for capy_dir in glob.glob(os.path.join(sys.argv[1], f'*capy*{mesh_to_plot}*')):
#     mesh = int(capy_dir.split('_')[-1])
#     results = np.genfromtxt(os.path.join(capy_dir, "results.csv"))
#     added_mass = results[:, 0]
#     damping = results[:, 1]

#     plt.plot(omega_range, added_mass, label=capy_dir)

# plt.legend()

plt.show()
