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

    omega = dataset.coords['omega']  # Range of frequencies in the dataset

    if 'excitation_force' not in dataset:
        dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']
    excitation = dataset['excitation_force'].sel(wave_direction=wave_direction)

    # SOLVE LINEAR SYSTEMS
    # Reorder dimensions of the arrays to be sure to solve the right system.
    H = H.transpose('omega', 'radiating_dof', 'influenced_dof')
    excitation = excitation.transpose('omega',  'influenced_dof')

    # Solve the linear systems (one for each value of omega)
    X = np.linalg.solve(H, excitation)

    return xr.DataArray(X, coords=[omega, dataset.coords['radiating_dof']], dims=['omega', 'radiating_dof'])

