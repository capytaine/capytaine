#!/usr/bin/env python
# coding: utf-8
"""Experimental function to compute the Response Amplitude Operator."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr
from capytaine.post_pro.impedance import rao_transfer_function

LOG = logging.getLogger(__name__)


def rao(dataset, wave_direction=0.0, dissipation=None, stiffness=None):
    """Response Amplitude Operator.

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'inertia_matrix' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    wave_direction: float, optional
        The direction of the incoming waves.
        Default: 0 radian (x direction)
    dissipation: array, optional
        An optional dissipation matrix (e.g. Power Take Off) to be included in the RAO.
        Default: none.
    stiffness: array, optional
        An optional stiffness matrix (e.g. mooring stiffness) to be included in the RAO.
        Default: none.

    Returns
    -------
    xarray DataArray
        The RAO as an array depending of omega and the degree of freedom.
    """

    # ASSEMBLE MATRICES
    H = rao_transfer_function(dataset, dissipation, stiffness)

    LOG.info("Compute RAO.")

    freq_dims = set(dataset.dims) & {'omega', 'period', 'wavelength', 'wavenumber'}
    if len(freq_dims) != 1:
        raise ValueError("The dataset provided to compute the RAO should one (and only one) dimension" +
                         " among the following: 'omega', 'period', 'wavelength' or 'wavenumber'\n" +
                         f"The received dataset has the following dimensions: {dataset.dims}")
    main_freq_type = freq_dims.pop()

    if 'excitation_force' not in dataset:
        dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']
    excitation = dataset['excitation_force'].sel(wave_direction=wave_direction)

    # SOLVE LINEAR SYSTEMS
    # Reorder dimensions of the arrays to be sure to solve the right system.
    H = H.transpose(main_freq_type, 'radiating_dof', 'influenced_dof')
    excitation = excitation.transpose(main_freq_type,  'influenced_dof')

    # Solve the linear systems (one for each value of omega)
    X = np.linalg.solve(H, excitation)

    return xr.DataArray(X, coords=[dataset.coords[main_freq_type], dataset.coords['radiating_dof']], dims=[main_freq_type, 'radiating_dof'])

