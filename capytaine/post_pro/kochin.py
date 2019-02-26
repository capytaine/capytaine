#!/usr/bin/env python
# coding: utf-8
"""
Computation of the Kochin function.
"""
# This file is part of "capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import numpy as np
import pandas as pd
import xarray as xr

from capytaine.bem.nemoh import Nemoh
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.io.xarray import _squeeze_dimensions, wavenumber_data_array


def compute_Kochin(result, theta, ref_point=(0.0, 0.0)):
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

    if result.sources is None:
        raise Exception(f"""The values of the sources of {result} cannot been found.
        They probably have not been stored by the solver because the option keep_details=True have not been set.
        Please re-run the resolution with this option.""")

    k = result.wavenumber
    h = result.depth

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


def kochin_data_array(results, theta_range, **kwargs):
    records = [dict(result.settings_dict, theta=theta, kochin=kochin)
               for result in results
               for theta, kochin in zip(theta_range, compute_Kochin(result, theta_range, **kwargs))]

    optional_vars = ['g', 'rho', 'body_name', 'water_depth']
    dimensions = ['omega', 'radiating_dof', 'theta']
    df = pd.DataFrame(records)
    df = df.set_index(optional_vars + dimensions)
    array = df.to_xarray()['kochin']

    array = _squeeze_dimensions(array, dimensions=optional_vars)

    return array


def kochin_dataset(floating_body, omega_range, theta_range):
    """Return a xarray dataset with all the Kochin functions for a given body and several frequencies."""
    solver = Nemoh()
    problems = [RadiationProblem(body=floating_body, radiating_dof=dof, omega=omega)
                for dof in floating_body.dofs for omega in omega_range]
    results = [solver.solve(problem, keep_details=True) for problem in sorted(problems)]

    kochin_dataset = kochin_data_array(results, theta_range)
    ds = xr.Dataset({'kochin': kochin_dataset})
    ds = ds.rename({'radiating_dof': 'dof'})

    ds.coords['nb_faces'] = ('mesh', [floating_body.mesh.nb_faces])

    # Add more stuffs to the dataset
    ds.coords['wavenumber'] = wavenumber_data_array(results)
    ds.coords['wavelength'] = 2*np.pi / ds.coords['wavenumber']

    return ds
