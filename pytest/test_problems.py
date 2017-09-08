#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np

from capytaine.problems import RadiationProblem
from capytaine.reference_bodies import DummyBody

def test_depth():
    dummy = DummyBody()
    assert RadiationProblem(dummy, free_surface=np.infty, sea_bottom=-np.infty).depth == np.infty
    assert RadiationProblem(dummy, free_surface=0.0, sea_bottom=-np.infty).depth == np.infty
    assert RadiationProblem(dummy, free_surface=0.0, sea_bottom=-1.0).depth == 1.0
    with pytest.raises(Exception):
        RadiationProblem(dummy, free_surface=0.0, sea_bottom=1.0)
