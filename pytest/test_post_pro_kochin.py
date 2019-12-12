#!/usr/bin/env python
# coding: utf-8
"""Test related to the computation of Kochin functions."""

import numpy as np
from numpy import pi
import xarray as xr

import capytaine as cpt


def test_fill_dataset_with_kochin_functions():
    sphere = cpt.Sphere(clip_free_surface=True)
    sphere.add_translation_dof(name="Heave")
    solver = cpt.BEMSolver()

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': [0.0, pi/2],
        'radiating_dof': ["Heave"],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin' in ds
    assert 'kochin_diffraction' not in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin'].sel(radiating_dof="Heave", theta=0.0),
                      ds['kochin'].sel(radiating_dof="Heave", theta=pi/2))

    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'radiating_dof': ["Heave"],
        'wave_direction': [-pi/2, 0.0],
        'theta': [0.0, pi/2],
    })
    ds = solver.fill_dataset(test_matrix, [sphere])
    assert 'theta' in ds.coords
    assert 'kochin' in ds
    assert 'kochin_diffraction' in ds

    # Because of the symmetries of the body
    assert np.isclose(ds['kochin_diffraction'].sel(wave_direction=-pi/2, theta=0.0),
                      ds['kochin_diffraction'].sel(wave_direction=0.0, theta=pi/2))
