"""This module implements objects describing a mesh on which the free surface elevation will be computed."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from itertools import product

import numpy as np

from capytaine.meshes.meshes import Mesh

LOG = logging.getLogger(__name__)


class FreeSurface():
    """A cartesian mesh on which the free surface elevation will be computed.

    Has a :code:`mesh` attribute to behave kind of like FloatingBody when
    building of the influence matrix.

    Parameters
    ----------
    x_range: Tuple[float, float], optional
        extreme values of the mesh in the x direction
    nx: int, optional
        number of cells in the x direction
    y_range: Tuple[float, float], optional
        extreme values of the mesh in the y direction
    ny: int, optional
        number of cells in the y direction
    name: string, optional
        a name for the free surface object


    .. todo:: Generalize to non-cartesian meshes.
              In particular, it could be of interest to build meshes having the
              same symmetry as a given floating body to speed up the
              construction of the influence matrix.

    .. seealso::

        :meth:`~capytaine.bem.nemoh.Nemoh.get_free_surface_elevation`
            The main function requiring a FreeSurface object.
    """
    def __init__(self, x_range=(-50.0, 50.0), nx=10, y_range=(-50.0, 50.0), ny=10, name=None):
        self.x_range = x_range
        self.nx = nx
        self.y_range = y_range
        self.ny = ny

        if name is None:
            self.name = f"free_surface_{next(Mesh._ids)}"
        else:
            self.name = name

        self.mesh = self._generate_mesh()

    def _generate_mesh(self):
        """Generate a 2D cartesian mesh."""
        nodes = np.zeros(((self.nx+1)*(self.ny+1), 3), dtype=float)
        panels = np.zeros((self.nx*self.ny, 4), dtype=int)

        X = np.linspace(*self.x_range, self.nx+1)
        Y = np.linspace(*self.y_range, self.ny+1)
        for i, (x, y, z) in enumerate(product(X, Y, [0.0])):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, self.nx), range(0, self.ny))):
            panels[k, :] = (j+i*(self.ny+1),
                            (j+1)+i*(self.ny+1),
                            (j+1)+(i+1)*(self.ny+1),
                            j+(i+1)*(self.ny+1))

        return Mesh(nodes, panels, name=f"{self.name}_mesh")

    @property
    def area(self):
        """The total area covered by the mesh."""
        return (np.abs(self.x_range[1] - self.x_range[0])
                * np.abs(self.y_range[1] - self.y_range[0]))

    def incoming_waves(self, problem: "DiffractionProblem") -> np.ndarray:
        """Free surface elevation of the undisturbed incoming waves
        for a given diffraction problem.
        Kept for legacy, but not recommended for use.
        """
        from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
        return airy_waves_free_surface_elevation(self, problem)
