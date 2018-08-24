#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

from capytaine.mesh.mesh import Mesh

import capytaine.mesh.mesh_clipper as mc

from capytaine.geometric_bodies import Sphere
from capytaine.tools.geometry import Plane


def test_clipper():
    mesh = Sphere().mesh.merge()

    plane = Plane()
    clipper = mc.MeshClipper(mesh, plane, assert_closed_boundaries=True)

    # for iter in range(50):
    #     thetax, thetay = np.random.rand(2)*2*np.pi
    #     plane.rotate(thetax, thetay)
    #     clipper.plane = plane
