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

from capytaine.problems import DiffractionProblem
from capytaine.tools.Airy_wave import Airy_wave_potential

LOG = logging.getLogger(__name__)


class FreeSurface():
    def __init__(self, x_range=(-50, 50), nx=10, y_range=(-50, 50), ny=10, name=None):
        self.x_range = x_range
        self.nx = nx
        self.y_range = y_range
        self.ny = ny

        self.nb_matrices_to_keep = 1

        if name is None:
            self.name = f"free_surface_{next(Mesh._ids)}"
        else:
            self.name = name

        self.mesh = self._generate_mesh()

    def _generate_mesh(self):
        nodes = np.zeros(((self.nx+1)*(self.ny+1), 3), dtype=np.float32)
        panels = np.zeros((self.nx*self.ny, 4), dtype=np.int)

        X = np.linspace(*self.x_range, self.nx+1)
        Y = np.linspace(*self.y_range, self.ny+1)
        for i, (x, y, z) in enumerate(product(X, Y, [0.0])):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, self.nx), range(0, self.ny))):
            panels[k, :] = (j+i*(self.ny+1), j+1+i*(self.ny+1), j+1+(i+1)*(self.ny+1), j+(i+1)*(self.ny+1))

        return Mesh(nodes, panels, name=f"{self.name}_mesh")

    @property
    def width(self):
        return np.abs(self.x_range[1] - self.x_range[0])

    @property
    def length(self):
        return np.abs(self.y_range[1] - self.y_range[0])

    @property
    def area(self):
        return self.width*self.length

    @staticmethod
    def with_same_symmetries_as(body):
        raise NotImplemented()

    # Tools for plotting of the free surface elevation

    def elevation_at_nodes(self, fs_faces: np.ndarray) -> np.ndarray:
        """From a free surface elevation computed at the center of the faces of the mesh,
        return a free surface elevation computed on the nodes of the mesh."""
        z_nodes = np.zeros((self.mesh.vertices.shape[0]), dtype=np.complex)
        faces_near_nodes = np.zeros((self.mesh.vertices.shape[0]), dtype=np.int)
        for i, vertices in enumerate(self.mesh.faces):
            for vertex in vertices:
                faces_near_nodes[vertex] += 1
                z_nodes[vertex] += fs_faces[i]
        z_nodes /= faces_near_nodes
        return z_nodes

    def incoming_waves(self, problem: DiffractionProblem) -> np.ndarray:
        """Free surface elevation of incoming wave for diffraction problem."""
        return 1j*problem.omega/problem.g * Airy_wave_potential(self.mesh.faces_centers, problem)
