"""Computation of the mean drift force."""
# Copyright (C) 2026 the Capytaine developers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import numpy as np
import xarray as xr

def far_field_mean_drift_force(X, dataset):
    """Compute the mean drift forces using far field formulation.
    Note that the forces are proportional to the square of the wave amplitude, but this implementation
    does not take into account the wave amplitude.

    Parameters
    ----------
    X : xarray DataArray
        The motion RAO.
    dataset : xarray Dataset
        This function supposes that variables named 'kochin_diffraction' and 'kochin_radiation' are in the dataset.

    Returns
    -------
    xarray Dataset
        The horizontal mean drift forces, depending on omega and the wave direction.
    """
    omega = dataset['omega']
    k = dataset['wavenumber']
    h = dataset['water_depth']
    beta = dataset['wave_direction']
    rho = dataset['rho']
    g = dataset['g']
    H_diff = dataset['kochin_diffraction']
    H_rad = dataset['kochin_radiation'] 
    theta_range = dataset['theta']

    if (theta_range.min() > 0) or (theta_range.max() < 2*np.pi):
        raise ValueError("Theta should takes values between 0 and 2pi")
    if np.any(beta < theta_range.min()) or np.any(beta > theta_range.max()):
        raise ValueError("The wave direction should be in the theta interval")
    if np.any(beta == theta_range.min()) or np.any(beta == theta_range.min()):
        raise ValueError("The wave direction should not be at the border of the theta interval, it is recommanded to extend the theta interval")

    H_rad_tot = sum(H_rad.sel(radiating_dof=d)*X.sel(radiating_dof=d) for d in X.radiating_dof)
    H_tot = np.exp(1j*np.pi/2)*(H_diff + H_rad_tot)
    H_beta = H_tot.interp(theta=beta)

    H_derivative = H_tot.differentiate("theta")
    H_derivative_beta = H_derivative.interp(theta=beta)

    coef1 = 2*np.pi*rho*omega
    if h == np.inf:
        coef2 = 2*np.pi*rho*k**2
    else:
        k0 = omega**2/g
        coef2 = 2*np.pi*rho*k*(k0*h)**2 / (h*((k*h)**2 - (k0*h)**2 + k0*h))

    dims = [d for d in X.dims if d != 'radiating_dof']
    coords = {c: X.coords[c] for c in X.coords if c != 'radiating_dof'}

    base = np.abs(H_tot.sel(theta=slice(0, 2*np.pi)))**2
    
    Fx = xr.DataArray(
        data=-coef1*np.cos(beta)*np.imag(H_beta) - coef2*(base * np.cos(base.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_surge'
        )

    Fy = xr.DataArray(
        data=-coef1*np.sin(beta)*np.imag(H_beta) - coef2*(base * np.sin(base.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_sway'
        )

    Mz = xr.DataArray(
        data=coef1/k*np.real(H_derivative_beta) - coef2/k*np.imag((np.conjugate(H_tot.sel(theta=slice(0, 2*np.pi))) * H_derivative.sel(theta=slice(0, 2*np.pi))).integrate("theta")),
        coords=coords,
        dims=dims,
        name='drift_force_yaw'
        )

    return xr.Dataset({Fx.name: Fx, Fy.name: Fy, Mz.name: Mz})


def _merge_far_field_mean_drift_variables(dataset):
    """Merge the three drift force components into a single variable with an influenced_dof dimension.

    Parameters
    ----------
    dataset : xarray Dataset
        The dataset returned by far_field_mean_drift_force, containing 'drift_force_surge',
        'drift_force_sway', and 'drift_force_yaw' variables.

    Returns
    -------
    xarray Dataset
        A dataset with a single 'drift_force' variable having an additional 'influenced_dof' dimension
        with coordinates ["Surge", "Sway", "Yaw"].
    """
    # Stack the three force components along a new dimension
    drift_force = xr.concat(
        [dataset['drift_force_surge'], dataset['drift_force_sway'], dataset['drift_force_yaw']],
        dim=xr.DataArray(['Surge', 'Sway', 'Yaw'], dims=['influenced_dof'], name='influenced_dof')
    )
    drift_force.name = 'drift_force'

    return xr.Dataset({'drift_force': drift_force})
