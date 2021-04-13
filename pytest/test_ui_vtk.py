#!/usr/bin/env python
# coding: utf-8

from warnings import warn

import capytaine as cpt

def test_animation_of_dofs():
    body = cpt.Sphere()
    body.add_translation_dof(name="Heave")
    try:
        animation = body.animate({"Heave": 0.2}, loop_duration=1.0)
        animation.embed_in_notebook()
    except ImportError as error:
        warn("VTK is not installed and thus has not been tested.")

