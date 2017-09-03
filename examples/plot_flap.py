#!/usr/bin/env python
# coding: utf-8

import glob
import numpy as np
import matplotlib.pyplot as plt

T_range, mu, nu = np.loadtxt("data/flap_mu_nu.tsv").T

plt.figure(1)
plt.plot(T_range, mu, linestyle="--", label="Reference added mass")

plt.figure(2)
plt.plot(T_range, nu, linestyle="--", label="Reference added damping")

for fichier in glob.glob("*.tsv"):
    mu, nu = np.loadtxt(fichier).T

    nb_cells = int(fichier.split('_')[1])

    plt.figure(1)
    plt.plot(
        T_range,
        mu,
        # color=f'{1-(i+1)/len(resolutions)}',
        label=f"Added mass ({nb_cells} cells)"
    )

    plt.figure(2)
    plt.plot(
        T_range,
        nu,
        # color=f'{1-(i+1)/len(resolutions)}',
        label=f"Added damping ({nb_cells} cells)"
    )

plt.legend()
plt.tight_layout()
plt.show()
