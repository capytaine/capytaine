#!/usr/bin/env python
# coding: utf-8
"""
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from functools import lru_cache

import numpy as np

from capytaine.bem.prony_decomposition import find_best_exponential_decomposition
import capytaine.bem.NemohCore as NemohCore

tabulated_integrals = lru_cache(maxsize=1)(NemohCore.initialize_green_wave.initialize_tabulated_integrals)
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

    def evaluate(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        Srankine, Vrankine = self.evaluate_rankine(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        if (free_surface == np.infty or
                (free_surface - sea_bottom == np.infty and wavenumber in (0, np.infty))):
            # No more terms in the Green function
            return Srankine, Vrankine

        Swave, Vwave = self.evaluate_wave(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        # The real valued matrices Srankine and Vrankine are automatically recasted as complex in the sum.
        Swave += Srankine
        Vwave += Vrankine
        return Swave, Vwave

    def evaluate_rankine(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        # RANKINE TERM

        S, V = NemohCore.green_rankine.build_matrices_rankine_source(
            mesh1.faces_centers, mesh1.faces_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
                                 )

        if free_surface == np.infty:
            # No free surface, no more terms in the Green function
            return S, V

        # REFLECTION TERM

        def reflect_vector(x):
            y = x.copy()
            y[:, 2] *= -1
            return y

        if free_surface - sea_bottom == np.infty:
            # INFINITE DEPTH
            def reflect_point(x):
                y = x.copy()
                # y[:, 2] = 2*free_surface - x[:, 2]
                y[:, 2] *= -1
                y[:, 2] += 2*free_surface
                return y
        else:
            # FINITE DEPTH
            def reflect_point(x):
                y = x.copy()
                # y[:, 2] = 2*sea_bottom - x[:, 2]
                y[:, 2] *= -1
                y[:, 2] += 2*sea_bottom
                return y

        Srefl, Vrefl = NemohCore.green_rankine.build_matrices_rankine_source(
            reflect_point(mesh1.faces_centers), reflect_vector(mesh1.faces_normals),
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
        )

        if free_surface - sea_bottom < np.infty or wavenumber == 0.0:
            S += Srefl
            V += Vrefl
        else:
            S -= Srefl
            V -= Vrefl

        return S, V

    def evaluate_wave(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
        depth = free_surface - sea_bottom
        if depth == np.infty:
            return NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, 0.0,
                *self.tabulated_integrals,
                np.empty(1), np.empty(1),  # Dummy arrays that won't actually be used by the fortran code.
                mesh1 is mesh2
            )
        else:
            a_exp, lamda_exp = find_best_exponential_decomposition(
                wavenumber*depth*np.tanh(wavenumber*depth),
                wavenumber*depth,
                method=self.finite_depth_prony_decomposition_method,
            )

            return NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, depth,
                *self.tabulated_integrals,
                lamda_exp, a_exp,
                mesh1 is mesh2
                )

