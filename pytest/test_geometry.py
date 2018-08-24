#!/usr/bin/env python
# coding: utf-8

from capytaine.tools.geometry import Plane
from capytaine.geometric_bodies import HorizontalCylinder

def test_plane():
    plane = Plane()
    plane.set_normal_from_angles(0., 0.)
    plane.get_normal_orientation_wrt_z()

    other_plane = Plane()
    assert other_plane == plane

    plane.flip_normal()
    try:
        plane.get_edge_intersection([-1, -1, -1], [-2, -2, -2])
    except RuntimeError:
        pass
    # plane.set_plane_parameters(1, 1, 1)
    plane.get_normal_orientation_wrt_z()
    plane.orthogonal_projection_on_plane(HorizontalCylinder().mesh.vertices)
    plane.get_origin()
