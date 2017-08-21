#!/usr/bin/env python
# coding: utf-8
"""
Solver for the BEM problem based on Nemoh's Green function.
"""

from itertools import product

import numpy as np
from numpy.linalg import solve

import NemohCore._Green as _Green


class Nemoh:
    """
    Solver for the BEM problem based on Nemoh's Green function.
    """
    def __init__(self):
        _Green.initialize_green_2.initialize_green()

    def build_matrices(self, bodies, wavenumber, omega, depth, g):

        if depth < np.infty:
            _Green.initialize_green_2.lisc(omega**2*depth/g, wavenumber*depth)

        # for body1, body2 in product(bodies, repeat=2):
        body = bodies[0]
        body1 = bodies[0]
        body2 = bodies[0]

        S = np.zeros((body.nb_faces, body.nb_faces), dtype=np.complex64)
        V = np.zeros((body.nb_faces, body.nb_faces), dtype=np.complex64)

        for i, j in product(range(body.nb_faces), range(body.nb_faces)):

            if depth == np.infty:
                SP1, SM1, VSP1, VSM1 = _Green.green_1.vav(
                    body1.faces_centers[i, :],
                    body2.vertices[body2.faces[j, :], :].T,
                    body2.faces_centers[j, :],
                    body2.faces_normals[j, :],
                    body2.faces_areas[j],
                    body2.faces_radiuses[j],
                    0.0,
                    -1)
                SP2, SM2, VSP2, VSM2 = _Green.green_2.vnsinfd(
                    wavenumber,
                    body1.faces_centers[i, :],
                    body2.faces_centers[j, :],
                    body2.faces_areas[j])

            else:
                SP1, SM1, VSP1, VSM1 = _Green.green_1.vav(
                    body1.faces_centers[i, :],
                    body2.vertices[body2.faces[j, :], :].T,
                    body2.faces_centers[j, :],
                    body2.faces_normals[j, :],
                    body2.faces_areas[j],
                    body2.faces_radiuses[j],
                    depth,
                    1)
                SP2, SM2, VSP2, VSM2 = _Green.green_2.vnsfd(
                    wavenumber,
                    body1.faces_centers[i, :],
                    body2.faces_centers[j, :],
                    body2.faces_areas[j],
                    depth)

            S[i, j] = SP1 + SP2
            V[i, j] = np.dot(body.faces_normals[i, :], VSP1 + VSP2)

        return S, V

    def solve(self, problem):

        S, V = self.build_matrices(
            problem.bodies,
            problem.wavenumber,
            problem.omega,
            problem.depth,
            problem.g
        )
        identity = np.identity(V.shape[0], dtype=np.float32)

        for body in problem.bodies:
            for dof in body.dof:
                sources = solve(V + identity/2, body.dof[dof])
                potential = S @ sources

                c = - problem.rho * potential @ \
                        (body.dof[dof] * body.faces_areas)
                added_mass = c.real
                added_damping = problem.omega * c.imag
        return added_mass, added_damping
