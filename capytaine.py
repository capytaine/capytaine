#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import matmul, solve

class RadiationProblem():
    def __init__(bodies, frequency=1.0, rho=1000.0, g=9.81):
        self.bodies = bodies
        self.rho = rho
        self.g = g
        self.frequency = frequency

    def solve(self):

        S, V = self.build_matrices()

        for body in self.bodies:
            for dof in body.dof:
                sources = solve(V, body.dof[dof])
                potential = matmul(S, problem.sources)
                print(potential)

    def build_matrices(self):
        pass

class Nemoh():


