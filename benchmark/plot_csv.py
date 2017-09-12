#!/usr/bin/env python
# coding: utf-8

import sys

import pandas as pd
import matplotlib.pyplot as plt

times = pd.DataFrame.from_csv(sys.argv[1])
ax = times.plot()
ax.set(xlabel='number of cells in mesh', ylabel='time (seconds)')
plt.grid()
plt.show()
