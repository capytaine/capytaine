#!/usr/bin/env python
# coding: utf-8
"""
Output of the result as BEMIO hdf5 file.

[1]: https://wec-sim.github.io/WEC-Sim/advanced_features.html#bemio

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import h5py
import numpy as np

nb_dofs_per_body = 6
nb_bodies = 1
nb_freq = 10
nb_wave_dir = 4


arr = np.random.rand(200)
with h5py.File('random.hdf5', 'w') as f:
    f.create_dataset('simulation_parameters/scaled', (1,))
    f.create_dataset('simulation_parameters/wave_dir', (nb_wave_dir,))
    f.create_dataset('simulation_parameters/water_depth', (1,))
    f.create_dataset('simulation_parameters/w', (nb_freq,))
    f.create_dataset('simulation_parameters/T', (nb_freq,))

    for i_body in range(nb_bodies):
        f.create_dataset(f"body{i_body}/properties/name", (1,), dtype='string')
        f.create_dataset(f"body{i_body}/properties/body_number", (1,), dtype=np.int)
        f.create_dataset(f"body{i_body}/properties/cg", (3,))  # Center of gravity
        f.create_dataset(f"body{i_body}/properties/cb", (3,))  # Center of buoyancy
        f.create_dataset(f"body{i_body}/properties/disp_vol", (1,))
        f.create_dataset(f"body{i_body}/hydro_coeffs/linear_restoring_stiffness", (nb_dofs_per_body, nb_dofs_per_body))
        f.create_dataset(f"body{i_body}/hydro_coeffs/excitation/re", (nb_freq, nb_wave_dir, nb_dofs_per_body))
        f.create_dataset(f"body{i_body}/hydro_coeffs/excitation/im", (nb_freq, nb_wave_dir, nb_dofs_per_body))
        f.create_dataset(f"body{i_body}/hydro_coeffs/added_mass/all", (nb_freq, nb_bodies*nb_dofs_per_body, nb_dofs_per_body))
        f.create_dataset(f"body{i_body}/hydro_coeffs/added_mass/inf_freq", (nb_bodies*nb_dofs_per_body, nb_dofs_per_body))
        f.create_dataset(f"body{i_body}/hydro_coeffs/radiation_damping/all", (nb_freq, nb_bodies*nb_dofs_per_body, nb_dofs_per_body))

with h5py.File('hinge_barge.nc', 'r') as f:
    print(list(f['added_mass'].shape))



