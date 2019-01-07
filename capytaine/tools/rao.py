#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


def rao(dataset, wave_angle=0.0, dissipation=None, stiffness=None):
    """Response Amplitude Operator.

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'mass' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    wave_angle: float, optional
        The direction of the incomming waves.
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
    LOG.info("Compute RAO.")

    # ASSEMBLE MATRICES
    omega = dataset.coords['omega']  # Range of frequencies in the dataset

    A = (-omega**2*(dataset['mass'] + dataset['added_mass'])
         + 1j*omega*dataset['radiation_damping']
         + dataset['hydrostatic_stiffness'])

    if dissipation is not None:
        A += 1j*omega*dissipation

    if stiffness is not None:
        A += stiffness

    excitation = dataset['excitation_force'].sel(angle=wave_angle)

    # SOLVE LINEAR SYSTEMS
    # Reorder dimensions of the arrays to be sure to solve the right system.
    A = A.transpose('omega', 'radiating_dof', 'influenced_dof')
    excitation = excitation.transpose('omega',  'influenced_dof')

    # Solve the linear systems (one for each value of omega)
    X = np.linalg.solve(A, excitation)

    return xr.DataArray(X, coords=[omega, dataset.coord['influenced_dof']], dims=['omega', 'influenced_dof'])

