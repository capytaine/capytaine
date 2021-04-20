#!/usr/bin/env python
# coding: utf-8
"""Experimental function to compute the impendance."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


def impedance(dataset, dissipation=None, stiffness=None):
    """Impedance (see, e.g., Falnes).
    
    @book{falnes2002ocean,
          title={Ocean Waves and Oscillating Systems: Linear Interactions Including Wave-Energy Extraction},
          author={Falnes, J.},
          isbn={9781139431934},
          url={https://books.google.com/books?id=bl1FyQjCklgC},
          year={2002},
          publisher={Cambridge University Press}
    }

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'mass' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    dissipation: array, optional
        An optional dissipation matrix (e.g. Power Take Off) to be included in the impedance.
        Default: none.
    stiffness: array, optional
        An optional stiffness matrix (e.g. mooring stiffness) to be included in the impedance.
        Default: none.

    Returns
    -------
    xarray DataArray
        The impedance as an array depending of omega and the degree of freedom.
    """

    if not hasattr(dataset, 'mass'):
        raise AttributeError('dataset must have a mass variable')

    if not hasattr(dataset, 'hydrostatic_stiffness'):
        raise AttributeError('dataset must have a hydrostatic_stiffness variable')

    LOG.info("Compute impedance.")

    # ASSEMBLE MATRICES
    omega = dataset.coords['omega']  # Range of frequencies in the dataset

    A = (-omega**2*(dataset['mass'] + dataset['added_mass'])
         + 1j*omega*dataset['radiation_damping']
         + dataset['hydrostatic_stiffness'])

    if dissipation is not None:
        A = A + 1j*omega*dissipation

    if stiffness is not None:
        A = A + stiffness

    return A

