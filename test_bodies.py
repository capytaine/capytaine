#!/usr/bin/env python
# coding: utf-8

import pytest
import numpy as np
from bodies import *

@pytest.mark.parametrize("size", np.linspace(1, 10, 2))
@pytest.mark.parametrize("ncells", [6 ,11, 16])
def test_parallelepiped(size, ncells):
    rp = RectangularParallelepiped(
            height=size, length=size, thickness=size, nh=ncells, nl=ncells, nth=ncells)
    assert np.allclose(
            rp.faces_areas, [(size/(ncells-1))**2] * rp.nb_faces, rtol=1e-3)
    assert np.allclose(
            rp.faces_radiuses, [size/(ncells-1)*np.sqrt(2)/2] * rp.nb_faces, rtol=1e-3)
