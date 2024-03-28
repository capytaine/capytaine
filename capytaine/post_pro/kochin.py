"""Computation of the Kochin function."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
import numpy as np

LOG = logging.getLogger(__name__)

def compute_kochin(result, theta, ref_point=(0.0, 0.0)):
    """Compute the far field coefficient

    Parameters
    ----------
    result: LinearPotentialFlowResult
        solved potential flow problem
    theta: float or 1-dim array of floats
        angles at which the coefficient is computed
    ref_point: couple of float, optional
        point of reference around which the far field coefficient is computed

    Returns
    -------
    H: same type as theta
        values of the Kochin function
    """

    if result.forward_speed != 0.0:
        LOG.warning("Kochin functions with forward speed have never been validated.")

    if result.sources is None:
        raise Exception(f"""The values of the sources of {result} cannot been found.
        They probably have not been stored by the solver because the option keep_details=True have not been set.
        Please re-run the resolution with this option.""")

    k = result.wavenumber
    h = result.water_depth

    # omega_bar.shape = (nb_faces, 2) @ (2, nb_theta)
    omega_bar = (result.body.mesh.faces_centers[:, 0:2] - ref_point) @ (np.cos(theta), np.sin(theta))

    if 0 <= k*h < 20:
        cih = np.cosh(k*(result.body.mesh.faces_centers[:, 2]+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*result.body.mesh.faces_centers[:, 2])

    # cih.shape = (nb_faces,)
    # omega_bar.T.shape = (nb_theta, nb_faces)
    # result.body.mesh.faces_areas.shape = (nb_faces,)
    zs = cih * np.exp(-1j * k * omega_bar.T) * result.body.mesh.faces_areas

    # zs.shape = (nb_theta, nb_faces)
    # result.sources.shape = (nb_faces,)
    return zs @ result.sources/(4*np.pi)
