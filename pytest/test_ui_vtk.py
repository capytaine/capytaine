#!/usr/bin/env python
# coding: utf-8

import capytaine as cpt

def test_animation_of_dofs():
    body = cpt.Sphere()
    body.add_translation_dof(name="Heave")
    animation = body.animate({"Heave": 0.2}, loop_duration=1.0)
    animation.embed_in_notebook()

