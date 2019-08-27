#!/usr/bin/env python
# coding: utf-8
"""
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from functools import lru_cache

import numpy as np

from capytaine.green_functions.prony_decomposition import find_best_exponential_decomposition
import capytaine.green_functions.Delhommeau_f90 as Delhommeau_f90

tabulated_integrals = lru_cache(maxsize=1)(Delhommeau_f90.initialize_green_wave.initialize_tabulated_integrals)
LOG = logging.getLogger(__name__)

class Delhommeau:
    """
    Parameters
    ----------
    tabulation_nb_integration_points: int, optional
        Number of points for the evaluation of the tabulated elementary integrals w.r.t. :math:`theta`
        used for the computation of the Green function (default: 251)
    finite_depth_prony_decomposition_method: string, optional
        The implementation of the Prony decomposition used to compute the finite depth Green function.

    Attributes
    ----------
    tabulated_integrals: 3-ple of arrays
        Tabulated integrals for the computation of the Green function.
    """
    def __init__(self,
                 tabulation_nb_integration_points=251,
                 finite_depth_prony_decomposition_method='fortran',
                 ):
        self.tabulated_integrals = tabulated_integrals(328, 46, tabulation_nb_integration_points)

        self.finite_depth_prony_decomposition_method = finite_depth_prony_decomposition_method

        self.exportable_settings = {
            'green_function': 'Delhommeau',
            'tabulation_nb_integration_points': tabulation_nb_integration_points,
            'finite_depth_prony_decomposition_method': finite_depth_prony_decomposition_method,
        }

        self._hash = hash(self.exportable_settings.values())

    def __hash__(self):
        return self._hash

    def evaluate(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):

        depth = free_surface - sea_bottom
        if free_surface == np.infty: # No free surface, only a single Rankine source term

            a_exp, lamda_exp = np.empty(1), np.empty(1)  # Dummy arrays that won't actually be used by the fortran code.

            coeffs = np.array((1.0, 0.0, 0.0))

        elif depth == np.infty:

            a_exp, lamda_exp = np.empty(1), np.empty(1)  # Idem

            if wavenumber == 0.0:
                coeffs = np.array((1.0, 1.0, 0.0))
            elif wavenumber == np.infty:
                coeffs = np.array((1.0, -1.0, 0.0))
            else:
                coeffs = np.array((1.0, -1.0, 1.0))

        else:  # Finite depth
            a_exp, lamda_exp = find_best_exponential_decomposition(
                wavenumber*depth*np.tanh(wavenumber*depth),
                wavenumber*depth,
                method=self.finite_depth_prony_decomposition_method,
            )
            if wavenumber == 0.0:
                raise NotImplementedError
            elif wavenumber == np.infty:
                raise NotImplementedError
            else:
                coeffs = np.array((1.0, 1.0, 1.0))

        # Main call to Fortran code
        return Delhommeau_f90.matrices.build_matrices(
            mesh1.faces_centers, mesh1.faces_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            wavenumber, 0.0 if depth == np.infty else depth,
            coeffs,
            *self.tabulated_integrals,
            lamda_exp, a_exp,
            mesh1 is mesh2
        )
