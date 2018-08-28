#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

from capytaine.geometric_bodies import Sphere


def test_clipper():
    mesh = Sphere(radius=5.0, ntheta=10).mesh.merge()
    aabb = mesh.axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, sea_bottom=-np.infty)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, sea_bottom=-1.0)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed

    # With CollectionOfMeshes (AxialSymmetry)
    mesh = Sphere(radius=5.0, ntheta=10).mesh
    aabb = mesh.merge().axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, sea_bottom=-np.infty)
    assert np.allclose(mesh.merge().axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, sea_bottom=-1.0)
    assert np.allclose(mesh.merge().axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed
