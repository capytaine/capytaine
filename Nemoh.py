#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import solve
from itertools import product
import NemohCore._Green as _Green

class Nemoh:
    def __init__(self):
        _Green.initialize_green_2.initialize_green()

    def build_matrices(self, bodies, wavenumber, omega, depth, g):

        if depth < np.infty:
            _Green.initialize_green_2.lisc(omega**2*depth/g, wavenumber*depth)

        # for body1, body2 in product(bodies, repeat=2):
        body=bodies[0]
        body1=bodies[0]
        body2=bodies[0]

        S = np.zeros((body.npanels, body.npanels), dtype=np.complex64)
        V = np.zeros((body.npanels, body.npanels), dtype=np.complex64)

        for i, j in product(range(body.npanels), range(body.npanels)):

            if depth == np.infty:
                SP1, SM1, VSP1, VSM1 = _Green.green_1.vav(
                        body1.center_of_mass[i, :],
                        body2.nodes[body2.panels[j, :], :].T,
                        body2.center_of_mass[j, :],
                        body2.normal[j, :],
                        body2.area[j],
                        body2.radius[j],
                        0.0,
                        -1)
                SP2, SM2, VSP2, VSM2 = _Green.green_2.vnsinfd(
                        wavenumber,
                        body1.center_of_mass[i, :],
                        body2.center_of_mass[j, :],
                        body2.area[j])

            else:
                SP1, SM1, VSP1, VSM1 = _Green.green_1.vav(
                        body1.center_of_mass[i, :],
                        body2.nodes[body2.panels[j, :], :].T,
                        body2.center_of_mass[j, :],
                        body2.normal[j, :],
                        body2.area[j],
                        body2.radius[j],
                        depth,
                        1)
                SP2, SM2, VSP2, VSM2 = _Green.green_2.vnsfd(
                        wavenumber,
                        body1.center_of_mass[i, :],
                        body2.center_of_mass[j, :],
                        body2.area[j],
                        depth)

            S[i, j] = SP1 + SP2
            V[i, j] = np.dot(body.normal[i, :], VSP1 + VSP2)

        return S, V

    def solve(self, problem):

        S, V = self.build_matrices(
            problem.bodies,
            problem.wavenumber,
            problem.omega,
            problem.depth,
            problem.g
        )
        Id = np.identity(V.shape[0], dtype=np.float32)

        for body in problem.bodies:
            for dof in body.dof:
                sources = solve(V + Id/2, body.dof[dof])
                potential = S @ sources
                
                c = - problem.rho * potential @ (body.dof[dof] * body.area)
                added_mass = c.real
                added_damping = problem.omega * c.imag
        return added_mass, added_damping

