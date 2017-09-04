#!/usr/bin/env python
# coding: utf-8

import pytest
from capytaine.problems import *

def test_depth():
    assert RadiationProblem([], free_surface=np.infty, sea_bottom=-np.infty).depth == np.infty
    assert RadiationProblem([], free_surface=0.0, sea_bottom=-np.infty).depth == np.infty
    assert RadiationProblem([], free_surface=0.0, sea_bottom=-1.0).depth == 1.0
    with pytest.raises(Exception):
        RadiationProblem([], free_surface=0.0, sea_bottom=1.0)
