#!/usr/bin/env python
# coding: utf-8
"""Draft of implementation for rigid bodies and multiple rigid bodies."""
# Copyright (C) 2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import numpy as np

from capytaine.bodies.bodies import FloatingBody


class RigidBody(FloatingBody):
    """A rigid body with exactly the six usual degrees of freedom."""

    def __init__(self, mesh=None, rotation_center=None, name=None):
        super().__init__(mesh=mesh, dofs=None, name=name)
        if rotation_center is not None:
            self.rotation_center = rotation_center
        self.add_all_rigid_body_dofs()

        self.mass = ...
        self.inertia_matrix = ...
        self.hydrostatic_stiffness = ...

    def _huygens_transport(self, point):
        p_g = self.rotation_center - point
        return self.mass * (np.dot(p_g, p_g) * np.eye(3) - np.outer(p_g, p_g))

    def join_bodies(*bodies):
        joined = super().__class__.join_bodies(*bodies)
        joined.mass = ...
        if all(isinstance(body, RigidBody) for body in bodies):
            joined.inertia_matrix = ...
            joined.hydrostatic_stiffness = ...

