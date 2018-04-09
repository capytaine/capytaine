#!/usr/bin/env python
# coding: utf-8
"""
Computation of the Kochin function.
"""

import numpy as np


def compute_Kochin(result, theta, ref_point=(0.0, 0.0)):
    """Compute the far field coefficient

    Parameters
    ----------
    result: LinearPotentialFlowResult
        solved potential flow problem
    theta: float
        angle at which the coefficient is computed
    ref_point: couple of float, optional
        point of reference around which the far field coefficient is computed

    Returns
    -------
    H: float
        value of the Kochin function
    """

    if result.sources is None:
        raise Exception(f"""The values of the sources of {result} cannot been found.
        They probably have not been stored by the solver because the option keep_details=True have not been set.
        Please re-run the resolution with this option.""")

    k = result.wavenumber
    h = result.depth

    omega_bar = (result.body.mesh.faces_centers[:, 0:2] - ref_point) @ (np.cos(theta), np.sin(theta))

    if 0 <= k*h < 20:
        cih = np.cosh(k*(result.body.mesh.faces_centers[:, 2]+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*result.body.mesh.faces_centers[:, 2])

    zs = cih[:] * np.exp(-1j * k * omega_bar[:]) * result.body.mesh.faces_areas[:]

    return result.sources @ zs/(4*np.pi)
