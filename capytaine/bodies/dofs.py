# Copyright (C) 2017-2022 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>


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
