#!/usr/bin/env python
# coding: utf-8
"""Tests for the computation of the Green function using Fortran routines."""

import pytest
from math import log10
import numpy as np
import capytaine as cpt
#----------------------------------------------------------------------------------#
gf = cpt.Delhommeau()

legacy_grid = gf.fortran_core.delhommeau_integrals.legacy_grid
scaled_nemoh3_grid = gf.fortran_core.delhommeau_integrals.scaled_nemoh3_grid
'''

for nr in [1000, 1500, 2000]:
    r_range = gf.fortran_core.delhommeau_integrals.default_r_spacing(nr, 100.0, 1)
    print(r_range[(r_range > 0.1) & (r_range < 1.0)])  # Print points between 0.1 and 1

'''
zmin = [-100,-251,-400]
@pytest.mark.parametrize("zmin", zmin)
def test_green_functions_z_spacing_zmin(zmin):
    """Test default_z_spacing"""

    nz = 124
    zRange = gf.fortran_core.delhommeau_integrals.default_z_spacing(nz, zmin, scaled_nemoh3_grid)

    dz = []
    for i in range(1,len(zRange)):
        dz.append(log10(abs(zRange[i]))-log10(abs(zRange[i-1])))
    dzComp = dz[-1]*np.ones_like(dz)

    assert abs(zmin-zRange[-1]) < 1e-5
    assert all([abs(a - b) < 1e-5 for a, b in zip(dz, dzComp)])


nz = [64,124,248]
@pytest.mark.parametrize("nz", nz)
def test_green_functions_z_spacing_nz(nz):
    """Test default_z_spacing"""
    zmin = -251
    zRange = gf.fortran_core.delhommeau_integrals.default_z_spacing(nz,zmin, scaled_nemoh3_grid)

    dz = []
    for i in range(1,len(zRange)):
        dz.append(log10(abs(zRange[i]))-log10(abs(zRange[i-1])))
    dzComp = dz[-1]*np.ones_like(dz)

    assert abs(zmin-zRange[-1]) < 1e-5
    assert all([abs(a - b) < 1e-5 for a, b in zip(dz, dzComp)])


nr = [100,328,480]
@pytest.mark.parametrize("nr", nr)
def test_green_functions_r_spacing_legacy_nr(nr):
    """Test default_r_spacing legacy model"""

    rmax = 4/3+abs(nr-32)/3
    rRange = gf.fortran_core.delhommeau_integrals.default_r_spacing(nr, rmax, legacy_grid)

    r_cal = []
    for i in range(1,nr+1):
        r_cal.append(min(10**((i-1)/5-6),4/3+abs(i-32)/3))

    assert all([abs(a - b) < 1e-4 for a, b in zip(rRange, r_cal)])


rmax = [50,100,200]
@pytest.mark.parametrize("rmax", rmax)
def test_green_functions_r_spacing_Nemoh_rmax(rmax):
    """Test default_r_spacing Nemoh model"""

    nr = 676
    c = 81
    r_logSpace = -np.log(10.0**(-10))/(c - 1.0)
    rRange = gf.fortran_core.delhommeau_integrals.default_r_spacing(nr, rmax, scaled_nemoh3_grid)

    r_cal = []
    for i in range(1,nr+1):
        if (i < c):
            r_cal.append(10**(-10)*np.exp(i*r_logSpace))
        else:
            r_cal.append((rmax-1)/(nr-c)*(i-c)+1)

    assert all([abs(a - b) < 1e-4 for a, b in zip(rRange, r_cal)])


nr = [676,1352]
@pytest.mark.parametrize("nr", nr)
def test_green_functions_r_spacing_Nemoh_nr(nr):
    """Test default_r_spacing Nemoh model"""

    nr0 = 676
    c = round(nr/nr0*81)
    rmax = 100
    r_logSpace = -np.log(10.0**(-10))/(c - 1.0)
    rRange = gf.fortran_core.delhommeau_integrals.default_r_spacing(nr, rmax, scaled_nemoh3_grid)

    r_cal = []
    for i in range(1,nr+1):
        if (i < c):
            r_cal.append(10**(-10)*np.exp(i*r_logSpace))
        else:
            r_cal.append((rmax-1)/(nr-c)*(i-c)+1)

    assert all([abs(a - b) < 1e-4 for a, b in zip(rRange, r_cal)])
