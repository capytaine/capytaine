#!/usr/bin/env python
# coding: utf-8

import pytest
import numpy as np

from capytaine.tools.geometry import Axis, Ox_axis, Oy_axis, Oz_axis, Plane, xOz_Plane, yOz_Plane, xOy_Plane
from capytaine.geometric_bodies import HorizontalCylinder


def test_axis():
    assert (2, 0, 0) in Ox_axis
    assert (0, 1, 0) not in Ox_axis


def test_plane():
    plane = Plane()
    other_plane = Plane()
    assert other_plane == plane


def test_plane_translations():
    translated_plane = xOz_Plane.translate(vector=(1, 0, 0), inplace=False)
    assert xOz_Plane is not translated_plane
    assert xOz_Plane == translated_plane

    translated_plane = xOz_Plane.translate(vector=(0, 1, 0), inplace=False)
    assert translated_plane.c == 1
    assert np.all(translated_plane.normal == xOz_Plane.normal)


def test_plane_rotations():
    rotated_plane = xOz_Plane.rotate(Oy_axis, angle=np.pi/12, inplace=False)
    assert xOz_Plane is not rotated_plane
    assert xOz_Plane == rotated_plane

    rotated_plane = xOz_Plane.rotate(Ox_axis, angle=np.pi/2, inplace=False)
    assert rotated_plane == xOy_Plane

    with pytest.raises(NotImplementedError):
        axis = Ox_axis.translate((0, 1, 0), inplace=False)
        xOz_Plane.rotate(axis, angle=1)

    with pytest.raises(NotImplementedError):
        plane = xOz_Plane.translate((0, 1, 0), inplace=False)
        plane.rotate(Oz_axis, angle=1)
