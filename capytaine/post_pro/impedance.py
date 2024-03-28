"""Computation of the impendance matrix."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

LOG = logging.getLogger(__name__)


def rao_transfer_function(dataset, dissipation=None, stiffness=None):
    """Complex-valued matrix used for the computation of the RAO.

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'inertia_matrix' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    dissipation: array, optional
        An optional dissipation matrix (e.g. Power Take Off) to be included in the transfer function.
        Default: none.
    stiffness: array, optional
        An optional stiffness matrix (e.g. mooring stiffness) to be included in the transfer function.
        Default: none.

    Returns
    -------
    xarray DataArray
        The matrix as an array depending of omega and the degrees of freedom.
    """

    if not hasattr(dataset, 'inertia_matrix'):
        raise AttributeError('Computing the impedance matrix requires an `inertia_matrix` matrix to be defined in the hydrodynamical dataset')

    if not hasattr(dataset, 'hydrostatic_stiffness'):
        raise AttributeError('Computing the impedance matrix requires an `hydrostatic_stiffness` matrix to be defined in the hydrodynamical dataset')

    if 'encounter_omega' in dataset.coords:
        omega = dataset.coords['encounter_omega']
    else:
        omega = dataset.coords['omega']

    # ASSEMBLE MATRICES
    H = (-omega**2*(dataset['inertia_matrix'] + dataset['added_mass'])
         - 1j*omega*dataset['radiation_damping']
         + dataset['hydrostatic_stiffness'])

    if dissipation is not None:
        H = H - 1j*omega*dissipation

    if stiffness is not None:
        H = H + stiffness

    return H


def impedance(dataset, dissipation=None, stiffness=None):
    """Complex-valued mechanical impedance matrix.
    See Falnes for more theoretical details::

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
        This function supposes that variables named 'inertia_matrix' and 'hydrostatic_stiffness' are in the dataset.
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
    if 'encounter_omega' in dataset.coords:
        omega = dataset.coords['encounter_omega']
    else:
        omega = dataset.coords['omega']
    return 1/(-1j * omega) * rao_transfer_function(dataset, dissipation, stiffness)
