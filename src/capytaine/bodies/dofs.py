# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

TRANSLATION_DOFS_DIRECTIONS = {"surge": (1, 0, 0), "sway": (0, 1, 0), "heave": (0, 0, 1)}
ROTATION_DOFS_AXIS = {"roll": (1, 0, 0), "pitch": (0, 1, 0), "yaw": (0, 0, 1)}


class RigidBodyDofsPlaceholder:
    """Pass an instance of this class to the FloatingBody initializer to initialize the 6 ridig body dofs."""

    def __init__(self, rotation_center=None):
        self.rotation_center = rotation_center

    def __str__(self):
        return "RigidBodyDofsPlaceholder()"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())


def rigid_body_dofs(rotation_center=None):
    return RigidBodyDofsPlaceholder(rotation_center=rotation_center)
