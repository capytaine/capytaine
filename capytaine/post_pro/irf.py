#!/usr/bin/env python
# coding: utf-8
"""Experimental function to compute the impendance."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr

LOG = logging.getLogger(__name__)


def rad_kernel(dataset, time_series):
	"""Impulse response function kernel (see, e.g., Falnes)
    
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
    time_series: numpy.ndarray
    	Time series on which to calculate the IRF kernel.

    Returns
    -------
    Kr: xarray DataArray
    	Impulse response function kernel
	"""

	Krf = lambda t: 2/np.pi*(np.cos(dataset.omega * t) * dataset.radiation_damping).integrate('omega')

	Kr1 = np.array(list(map(Krf, time_series)))
	Kr = xr.DataArray(Kr1, dims=['time','radiating_dof','influenced_dof'],
	                  name='impulse_response_function',
	                  attrs=dataset.attrs,
	                  coords=dict(time=time_series,
	                                   radiating_dof=dataset.radiating_dof,
	                                   influenced_dof=dataset.influenced_dof
	                                   ))
	return Kr
