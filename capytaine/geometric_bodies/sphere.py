#!/usr/bin/env python
# coding: utf-8
"""
Generate meshed free surface.

This file is part of "capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging
from itertools import product

import numpy as np

from meshmagick.mesh import Mesh
from capytaine.bodies import FloatingBody
# from capytaine.symmetries import xOz_Plane, yOz_Plane, ReflectionSymmetry, TranslationalSymmetry, AxialSymmetry

LOG = logging.getLogger(__name__)


class Sphere(FloatingBody):
    pass
