#!/usr/bin/env python
# coding: utf-8

import numpy as np


def test_MaxLengthDict():
    from capytaine.tools.max_length_dict import MaxLengthDict

    dc = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=4)
    assert dc.__max_length__ == 4
    assert list(dc.keys()) == ['a', 'b', 'c']
    dc['d'] = 8
    assert list(dc.keys()) == ['a', 'b', 'c', 'd']
    dc['e'] = 2 # drop 'a'
    assert list(dc.keys()) == ['b', 'c', 'd', 'e']
    dc['b'] = 7 # move 'b' to front
    assert list(dc.keys()) == ['c', 'd', 'e', 'b']
    dc['f'] = 4 # drop 'c'
    assert list(dc.keys()) == ['d', 'e', 'b', 'f']
    dc.update({'g':6, 'h':9})
    assert list(dc.keys()) == ['b', 'f', 'g', 'h']

    dc2 = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=1)
    assert dc2.__max_length__ == 1
    assert list(dc2.keys()) == ['c']

    dc3 = MaxLengthDict({'a':1, 'b':5, 'c':3}, max_length=0)
    assert dc3 == {}
    dc3['d'] = 8
    assert dc3 == {}


def test_Froude_Krylov():
    from capytaine.tools.Airy_wave import Froude_Krylov_force
    from capytaine.reference_bodies import generate_clever_sphere
    from capytaine.problems import DiffractionProblem

    sphere = generate_clever_sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)

    problem = DiffractionProblem(body=sphere, omega=1.0, sea_bottom=-np.infty)
    assert np.isclose(Froude_Krylov_force(problem)['Heave'], 27596, rtol=1e-3)

    problem = DiffractionProblem(body=sphere, omega=2.0, sea_bottom=-np.infty)
    assert np.isclose(Froude_Krylov_force(problem)['Heave'], 22491, rtol=1e-3)

    problem = DiffractionProblem(body=sphere, omega=1.0, sea_bottom=-10.0)
    assert np.isclose(Froude_Krylov_force(problem)['Heave'], 27610, rtol=1e-3)


def test_Airy():
    """Compare finite depth Airy wave expression with results from analytical
    expression"""
    from capytaine.problems import DiffractionProblem
    from itertools import product
    from capytaine.tools.Airy_wave import Airy_wave_velocity, Airy_wave_potential

    try:
        import sympy as sp
        from sympy.abc import t
        from sympy.physics.vector import ReferenceFrame
        from sympy.physics.vector import gradient, divergence

        R = ReferenceFrame('R')
        x, y, z = R[0], R[1], R[2]
        Phi, k, h, g, rho = sp.symbols("Phi, k, h, g, rho")

        omega = sp.sqrt(g*k*sp.tanh(k*h))
        phi = g/omega * sp.cosh(k*(z+h))/sp.cosh(k*h) * sp.sin(k*x - omega*t)
        u = gradient(phi, R)
        p = -rho*Phi.diff(t)

        for depth in np.linspace(100.0, 10.0, 2):
            for omega in np.linspace(0.5, 4.0, 2):

                dp = DiffractionProblem(free_surface=0.0, sea_bottom=-depth, omega=omega)

                for t_val, x_val, y_val, z_val in product(np.linspace(0.0, 1.0, 2),
                                                          np.linspace(-10, 10, 3),
                                                          np.linspace(-10, 10, 3),
                                                          np.linspace(-10, 0, 3)):

                    parameters = {t: t_val, x: x_val, y: y_val, z: z_val,
                                  omega: dp.omega, k: dp.wavenumber,
                                  h: dp.depth, g: dp.g, rho: dp.rho}

                    phi_num = Airy_wave_potential(np.array((x_val, y_val, z_val)), dp)
                    assert np.isclose(float(phi.subs(parameters)),
                                      np.real(phi_num*np.exp(-1j * dp.omega * t_val)),
                                      rtol=1e-3)

                    u_num = Airy_wave_velocity(np.array((x_val, y_val, z_val)), dp)
                    assert np.isclose(float(u.dot(R.x).subs(parameters)),
                                      np.real(u_num[0]*np.exp(-1j * dp.omega * t_val)),
                                      rtol=1e-3)
                    assert np.isclose(float(u.dot(R.y).subs(parameters)),
                                      np.real(u_num[1]*np.exp(-1j * dp.omega * t_val)),
                                      rtol=1e-3)
                    assert np.isclose(float(u.dot(R.z).subs(parameters)),
                                      np.real(u_num[2]*np.exp(-1j * dp.omega * t_val)),
                                      rtol=1e-3)

    except ImportError:
        print("Airy wave not tested with sympy.")
