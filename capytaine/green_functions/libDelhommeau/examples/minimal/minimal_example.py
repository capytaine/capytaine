#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import capytaine as cpt


def print_matrix(x):
    for i in range(2):
        print(" %+19.12E %+19.12E %+19.12E %+19.12E" %
              (np.real(x[i, 0]), np.imag(x[i, 0]),
               np.real(x[i, 1]), np.imag(x[i, 1])))
    return


if __name__ == '__main__':

    print(' -- Run libdelhommeau/examples/minimal/minimal_example.py')

    vert = np.array([0.0, 0.0, -1.0,
                     1.0, 0.0, -1.0,
                     1.0, 1.0, -1.0,
                     0.0, 1.0, -1.0,
                     1.0, 0.0, -1.0,
                     2.0, 0.0, -1.0,
                     2.0, 1.0, -1.0,
                     1.0, 1.0, -1.0]).reshape((8, 3))

    faces = np.arange(0, 8).reshape(2, 4)
    mesh = cpt.Mesh(vert, faces)

    S, K = cpt.Delhommeau().evaluate(mesh, mesh, water_depth=np.inf, free_surface=np.inf)
    print(" Rankine part: S")
    print_matrix(S)
    print(" Rankine part: K")
    print_matrix(K)

    print(" k=1.0, h=inf: S")
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, water_depth=np.inf, wavenumber=1.0)
    print_matrix(S)
    print(" k=1.0, h=inf: K")
    print_matrix(K)

    S, K = cpt.Delhommeau().evaluate(mesh, mesh, water_depth=np.inf, wavenumber=2.0)
    print(" k=2.0, h=inf: S")
    print_matrix(S)
    print(" k=2.0, h=inf: K")
    print_matrix(K)

    print(" k=1.0, h=2.0: S")
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, water_depth=2., wavenumber=1.0)
    print_matrix(S)
    print(" k=1.0, h=2.0: K")
    print_matrix(K)

    S, K = cpt.Delhommeau().evaluate(mesh, mesh, water_depth=2., wavenumber=2.0)
    print(" k=2.0, h=2.0: S")
    print_matrix(S)
    print(" k=2.0, h=2.0: K")
    print_matrix(K)
