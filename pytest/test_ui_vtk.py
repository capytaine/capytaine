#!/usr/bin/env python
# coding: utf-8

import pytest
import capytaine as cpt

try:
    import vtk
except ImportError:
    vtk = None

@pytest.mark.skipif(vtk is None,
                    reason='vtk is not installed')
def test_animation_of_dofs():
    body = cpt.Sphere()
    body.add_translation_dof(name="Heave")
    animation = body.animate({"Heave": 0.2}, loop_duration=1.0)
    animation.embed_in_notebook()
