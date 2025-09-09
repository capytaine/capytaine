"""Experimental function to compute the Response Amplitude Operator."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr
from capytaine.post_pro.impedance import rao_transfer_function

LOG = logging.getLogger(__name__)


def rao(dataset, wave_direction=None, dissipation=None, stiffness=None):
    """Response Amplitude Operator.

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'inertia_matrix' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    wave_direction: float, optional
        Select a wave directions for the computation. (Not recommended, kept for legacy.)
        Default: all wave directions in the dataset.
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
    fex = dataset.excitation_force

    LOG.info("Compute RAO.")

    # SOLVE LINEAR SYSTEMS
    # Match dimensions of the arrays to be sure to solve the right systems.
    H, fex = xr.broadcast(H, fex, exclude=["radiating_dof", "influenced_dof"])
    H = H.transpose(..., 'radiating_dof', 'influenced_dof')
    fex = fex.transpose(...,  'influenced_dof')

    if wave_direction is not None:  # Legacy behavior for backward compatibility
        H = H.sel(wave_direction=wave_direction)
        fex = fex.sel(wave_direction=wave_direction)

    # Solve and add coordinates
    rao_dims = [d for d in H.dims if d != 'influenced_dof']
    rao_coords = {c: H.coords[c] for c in H.coords if c != 'influenced_dof'}
    rao = xr.DataArray(np.linalg.solve(H.values, fex.values[..., np.newaxis])[..., 0], coords=rao_coords, dims=rao_dims)

    return rao
