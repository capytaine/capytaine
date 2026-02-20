"""Computation of the mean drift force."""
# Copyright (C) 2026 the Capytaine developpers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import numpy as np
import xarray as xr

def far_field_formulation(X, dataset):
    omega = dataset['omega']
    k = dataset['wavenumber']
    h = dataset['water_depth']
    beta = dataset['wave_direction']
    rho = dataset['rho']
    g = dataset['g']
    H_diff = dataset['kochin_diffraction']
    H_rad = dataset['kochin_radiation'] 

    H_rad_tot = sum(H_rad.sel(radiating_dof=d)*X.sel(radiating_dof=d) for d in X.radiating_dof) 
    H_tot = np.exp(1j*np.pi/2)*(H_diff + H_rad_tot)
    H_beta = H_tot.sel(theta=beta) #assumes that wave direction is in the list of angle passed to the Kochin function 

    coef1 = 2*np.pi*rho*omega
    if h == np.inf:
        coef2 = 2*np.pi*rho*k**2
    else:
        k0 = omega**2/g 
        coef2 = 2*np.pi*rho*k*(k0*h)**2 / (h*((k*h)**2 - (k0*h)**2 + k0*h))

    dims = [d for d in X.dims if d != 'radiating_dof']
    coords = {c: X.coords[c] for c in X.coords if c != 'radiating_dof'}

    base = np.abs(H_tot)**2
    
    Fx = xr.DataArray(
        data=-coef1*np.cos(beta)*np.imag(H_beta) - coef2*(base * np.cos(H_tot.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_surge'
        )
    
    Fy = xr.DataArray(
        data=-coef1*np.sin(beta)*np.imag(H_beta) - coef2*(base * np.sin(H_tot.theta)).integrate("theta"),
        coords=coords,
        dims=dims,
        name='drift_force_sway'
        )
    
    return xr.Dataset({Fx.name: Fx, Fy.name: Fy})