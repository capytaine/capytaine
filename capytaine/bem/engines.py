#!/usr/bin/env python
# coding: utf-8
"""
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

from functools import lru_cache

import numpy as np

from capytaine.matrices import linear_solvers
from capytaine.matrices.builders import identity_like

LOG = logging.getLogger(__name__)

class BasicEngine:
    """
    Parameters
    ----------
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver
        (available: "direct" and "gmres", the latter is the default choice)
        or by passing directly a solver function.
    """
    available_linear_solvers = {'direct': linear_solvers.solve_directly,
                                'gmres': linear_solvers.solve_gmres}

    def __init__(self,
                 linear_solver='gmres',
                 matrix_cache_size=1,
                 ):

        if linear_solver in self.available_linear_solvers:
            self.linear_solver = self.available_linear_solvers[linear_solver]
        else:
            self.linear_solver = linear_solver

        if matrix_cache_size > 0:
            self.build_matrices = lru_cache(maxsize=matrix_cache_size)(self.build_matrices)

        self.exportable_settings = {
            'engine': 'BasicEngine',
            'matrix_cache_size': matrix_cache_size,
            'linear_solver': str(linear_solver),
        }

    def build_matrices(self,
                       mesh1, mesh2, free_surface, sea_bottom, wavenumber,
                       green_function):
        """ """

        S, V = green_function.evaluate(
            mesh1, mesh2, free_surface, sea_bottom, wavenumber,
        )

        return S, V + identity_like(V)/2

    def build_S_matrix_for_reconstruction(self, problem, mesh, green_function):
            S, _ = green_function.evaluate(
                mesh,
                problem.body.mesh,
                free_surface=problem.free_surface,
                sea_bottom=problem.sea_bottom,
                wavenumber=problem.wavenumber
            )
            return S


from capytaine.bem.hierarchical_toeplitz_matrices import hierarchical_toeplitz_matrices

class HierarchicalToeplitzMatrices:
    """

    Parameters
    ----------
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """
    def __init__(self,
                 ACA_distance=np.infty,
                 ACA_tol=1e-2,
                 matrix_cache_size=1,
                 ):

        if matrix_cache_size > 0:
            self.build_matrices = lru_cache(maxsize=matrix_cache_size)(self.build_matrices)

        self.ACA_distance = ACA_distance
        self.ACA_tol = ACA_tol

        self.linear_solver = linear_solvers.solve_gmres

        self.exportable_settings = {
            'engine': 'HierarchicalToeplitzMatrices',
            'ACA_distance': ACA_distance,
            'ACA_tol': ACA_tol,
            'matrix_cache_size': matrix_cache_size,
        }

    def build_matrices(self,
                       mesh1, mesh2, free_surface, sea_bottom, wavenumber,
                       green_function):
        """ """
        S, V = hierarchical_toeplitz_matrices(
            green_function.evaluate,
            ACA_tol=self.ACA_tol,
            ACA_distance=self.ACA_distance,
            dtype=np.complex128
        )(
            mesh1, mesh2, free_surface, sea_bottom, wavenumber,
        )

        return S, V + identity_like(V)/2

