#!/usr/bin/env python
# coding: utf-8
"""Output of the result as BEMIO hdf5 file.

[1]: https://wec-sim.github.io/WEC-Sim/advanced_features.html#bemio
"""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import h5py
import numpy as np


def to_bemio_file(dataset, body, filepath):
    """Write a dataset to BEMIO hdf5 file format.
    For the moment only a single body with 6 rigid body dofs is supported.

    Parameters
    ----------
    dataset: xarray.Dataset
        the dataset returned by capytaine.results.assemble_dataset
    body: FloatingBody
        the studied floating body
    filepath:
        path of the hdf5 file
    """
    nb_dofs_per_body = 6
    nb_bodies = 1
    nb_freq = len(dataset.coords['omega'])
    nb_wave_dir = len(dataset.coords['wave_direction']) if 'wave_direction' in dataset.coords else 0

    with h5py.File(filepath, 'w') as h5file:
        # h5file.create_dataset('simulation_parameters/scaled', (1,))

        h5file.create_dataset('simulation_parameters/wave_dir', (nb_wave_dir,))
        if nb_wave_dir > 0:
            h5file['simulation_parameters/wave_dir'][:] = dataset.coords['wave_direction'].values

        h5file.create_dataset('simulation_parameters/water_depth', (1,))
        h5file['simulation_parameters/water_depth'][0] = dataset.attrs['depth']

        h5file.create_dataset('simulation_parameters/w', (nb_freq,))
        h5file['simulation_parameters/w'][:] = dataset.coords['omega'].values

        h5file.create_dataset('simulation_parameters/T', (nb_freq,))
        h5file['simulation_parameters/T'][:] = 2*np.pi/dataset.coords['omega'].values

        for i_body in range(1, nb_bodies+1):
            # h5file.create_dataset(f"body{i_body}/properties/name", (len(body.name),), dtype='S10')
            # h5file[f"body{i_body}/properties/name"] = bytes(body.name)

            h5file.create_dataset(f"body{i_body}/properties/body_number", (1,), dtype=np.int)
            h5file[f"body{i_body}/properties/body_number"][0] = i_body

            h5file.create_dataset(f"body{i_body}/properties/cg", (3,))  # Center of gravity
            h5file[f"body{i_body}/properties/cg"][:] = body.center_of_gravity

            h5file.create_dataset(f"body{i_body}/properties/cb", (3,))  # Center of buoyancy
            h5file[f"body{i_body}/properties/cb"][:] = body.center_of_buoyancy

            h5file.create_dataset(f"body{i_body}/properties/disp_vol", (1,))
            h5file[f"body{i_body}/properties/disp_vol"][:] = body.displacement_volume

            h5file.create_dataset(f"body{i_body}/hydro_coeffs/linear_restoring_stiffness", (nb_dofs_per_body, nb_dofs_per_body))
            # h5file[f"body{i_body}/hydro_coeffs/linear_restoring_stiffness"][:] = body.

            h5file.create_dataset(f"body{i_body}/hydro_coeffs/excitation/re", (nb_freq, nb_wave_dir, nb_dofs_per_body))
            h5file.create_dataset(f"body{i_body}/hydro_coeffs/excitation/im", (nb_freq, nb_wave_dir, nb_dofs_per_body))
            if nb_wave_dir > 0:
                dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']
                h5file[f"body{i_body}/hydro_coeffs/excitation/re"][:] = np.real(dataset['excitation_force'].values)
                h5file[f"body{i_body}/hydro_coeffs/excitation/im"][:] = np.imag(dataset['excitation_force'].values)

            h5file.create_dataset(f"body{i_body}/hydro_coeffs/added_mass/all", (nb_freq, nb_bodies*nb_dofs_per_body, nb_dofs_per_body))
            h5file[f"body{i_body}/hydro_coeffs/added_mass/all"][:] = dataset['added_mass'].values

            h5file.create_dataset(f"body{i_body}/hydro_coeffs/added_mass/inf_freq", (nb_bodies*nb_dofs_per_body, nb_dofs_per_body))
            h5file[f"body{i_body}/hydro_coeffs/added_mass/inf_freq"][:] = dataset['added_mass'].sel(omega=max(dataset.coords['omega'].values)).values

            h5file.create_dataset(f"body{i_body}/hydro_coeffs/radiation_damping/all", (nb_freq, nb_bodies*nb_dofs_per_body, nb_dofs_per_body))
            h5file[f"body{i_body}/hydro_coeffs/radiation_damping/all"][:] = dataset['radiation_damping'].values


