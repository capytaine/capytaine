#!/usr/bin/env python
# coding: utf-8
"""
Generate meshed free surface.

This file is part of "capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
from itertools import product

from attr import attrs, attrib

import numpy as np

from meshmagick.mesh import Mesh
from capytaine.bodies import FloatingBody

LOG = logging.getLogger(__name__)


# class FreeSurface():
#     width = attrib()
#     length = attrib()
#     nw = attrib()
#     nl = attrib()

#     @staticmethod
#     def with_same_symmetries_as(body):
#         pass


def generate_free_surface(width=100, length=100, nw=10, nl=10, name=None):
    """Special body for the meshing of the free surface."""
    X = np.linspace(-width/2, width/2, nw+1)
    Y = np.linspace(-length/2, length/2, nl+1)

    nodes = np.zeros(((nw+1)*(nl+1), 3), dtype=np.float32)
    panels = np.zeros((nw*nl, 4), dtype=np.int)

    for i, (x, y, z) in enumerate(product(X, Y, [0.0])):
        nodes[i, :] = x, y, z

    for k, (i, j) in enumerate(product(range(0, nw), range(0, nl))):
        panels[k, :] = (j+i*(nl+1), j+1+i*(nl+1), j+1+(i+1)*(nl+1), j+(i+1)*(nl+1))

    if name is None:
        name = f"free_surface_{next(Mesh._ids)}"
    return FloatingBody(Mesh(nodes, panels), name=name)
