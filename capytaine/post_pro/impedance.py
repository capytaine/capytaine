#!/usr/bin/env python
# coding: utf-8
"""Computation of the impendance matrix."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


def impedance(dataset, dissipation=None, stiffness=None):
    """Complex-valued mechanical impedance matrix

    See Falnes for more theoretical details.
    Note that, unlike Falnes, we define here the impedance with respect to the position, that is `force = impedance x position`.
    For the impedance with respect to the velocity, see the function :code:`velocity_impedance` in the same module.

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
        The impedance as an array depending of omega and the degrees of freedom.
    """

    if not hasattr(dataset, 'mass'):
        raise AttributeError('Computing the impedance matrix requires a :code:`mass` matrix to be defined in the hydrodynamical dataset')

    if not hasattr(dataset, 'hydrostatic_stiffness'):
        raise AttributeError('Computing the impedance matrix requires a :code:`hydrostatic_stiffness` matrix to be defined in the hydrodynamical dataset')

    LOG.info("Compute impedance.")

    # ASSEMBLE MATRICES
    omega = dataset.coords['omega']  # Range of frequencies in the dataset

    A = (-omega**2*(dataset['mass'] + dataset['added_mass'])
         - 1j*omega*dataset['radiation_damping']
         + dataset['hydrostatic_stiffness'])

    if dissipation is not None:
        A = A - 1j*omega*dissipation

    if stiffness is not None:
        A = A + stiffness

    return A


def position_impedance(dataset, dissipation=None, stiffness=None):
    """Alias for the :code:`impedance` function."""
    return impedance(dataset, dissipation, stiffness)


def velocity_impedance(dataset, dissipation=None, stiffness=None):
    """Impedance with respect to the velocity, that is `force = velocity_impedance x velocity`.
    Takes the same arguments as :code:`impedance` which returns the impedance with respect to the position.
    """
    return 1/(-1j * dataset.coords["omega"]) * impedance(dataset, dissipation, stiffness)
