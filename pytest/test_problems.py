#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.problems import *
from capytaine.reference_bodies import DummyBody

dummy = DummyBody()

def test_depth():
    assert PotentialFlowProblem(dummy, free_surface=np.infty, sea_bottom=-np.infty).depth == np.infty
    assert PotentialFlowProblem(dummy, free_surface=0.0, sea_bottom=-np.infty).depth == np.infty
    assert PotentialFlowProblem(dummy, free_surface=0.0, sea_bottom=-1.0).depth == 1.0

    with pytest.raises(Exception):
        PotentialFlowProblem(dummy, free_surface=0.0, sea_bottom=1.0)

# def test_Airy():
#     dp = DiffractionProblem(dummy, free_surface=0.0, sea_bottom=-np.infty, omega=1.0)
#     print(dp.airy_wave((0.0, 0.0, -5.0)))

