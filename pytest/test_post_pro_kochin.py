#!/usr/bin/env python
# coding: utf-8
"""Test related to the computation of Kochin functions."""

import numpy as np
from numpy import pi
import xarray as xr

import capytaine as cpt

def test_kochin_sphere():
    sphere = cpt.Sphere(clip_free_surface=True)
    sphere.add_translation_dof(name="Heave")

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'wave_direction': [0.0, pi/2],
        'radiating_dof': ["Heave"],
        'theta': [0.0, pi/2, pi, 3*pi/2],
    })
    ds = cpt.BEMSolver().fill_dataset(test_matrix, [sphere])

    assert 'kochin' in ds
    assert 'kochin_diffraction' in ds
    # return ds

# ds = test_kochin_sphere()
# print(ds)
