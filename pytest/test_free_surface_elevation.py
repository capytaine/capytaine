#!/usr/bin/env python
# coding: utf-8
"""Tests related to post-processing computation of potential and free surface elevation."""

import pytest
import numpy as np
import capytaine as cpt


def test_potential():
    sphere = cpt.Sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_translation_dof(name="Heave")
    solver = cpt.BEMSolver()
    result = solver.solve(cpt.RadiationProblem(body=sphere, omega=1.0), keep_details=True)
    free_surface = cpt.FreeSurface(x_range=(-100, 100), nx=5, y_range=(-100, 100), ny=5)
    eta = solver.get_potential_on_mesh(result, free_surface.mesh, chunk_size=3)

