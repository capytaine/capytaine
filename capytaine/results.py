#!/usr/bin/env python
# coding: utf-8
"""
Definition of the objects storing the results of a LinearPotentialFlowProblem.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import logging

from attr import attrs, attrib, Factory
import numpy as np

from capytaine.tools.Airy_wave import Froude_Krylov_force


LOG = logging.getLogger(__name__)


@attrs(slots=True)
class LinearPotentialFlowResult:
    problem = attrib()

    sources = attrib(default=None, init=False, repr=False)
    potential = attrib(default=None, init=False, repr=False)

    def __getattr__(self, name):
        try:
            return getattr(self.problem, name)
        except AttributeError:
            raise AttributeError(f"{self.__class__} does not have a attribute named {name}.")


@attrs(slots=True)
class DiffractionResult(LinearPotentialFlowResult):
    forces = attrib(default=Factory(dict), init=False, repr=False)

    def __str__(self):
        parameters = [f"body={self.body.name}, omega={self.omega:.3f}, depth={self.depth}, angle={self.angle}, "]
        if not self.free_surface == 0.0:
            parameters.append(f"free_surface={self.free_surface}, ")
        if not self.g == 9.81:
            parameters.append(f"g={self.g}, ")
        if not self.rho == 1000:
            parameters.append(f"rho={self.rho}, ")
        return "DiffractionResult(" + ''.join(parameters)[:-2] + ")"

    def store_force(self, dof, force):
        self.forces[dof] = 1j*self.omega*force

    def records(self):
        FK = Froude_Krylov_force(self.problem)
        return [dict(omega=self.omega, influenced_dof=dof, angle=self.angle,
                     diffraction_force=self.forces[dof], Froude_Krylov_force=FK[dof])
                for dof in self.influenced_dofs]


@attrs(slots=True)
class RadiationResult(LinearPotentialFlowResult):
    added_masses = attrib(default=Factory(dict), init=False, repr=False)
    radiation_dampings = attrib(default=Factory(dict), init=False, repr=False)

    def __str__(self):
        parameters = [f"body={self.body.name}, omega={self.omega:.3f}, depth={self.depth}, radiating_dof={self.radiating_dof}, "]
        if not self.free_surface == 0.0:
            parameters.append(f"free_surface={self.free_surface}, ")
        if not self.g == 9.81:
            parameters.append(f"g={self.g}, ")
        if not self.rho == 1000:
            parameters.append(f"rho={self.rho}, ")
        return "RadiationResult(" + ''.join(parameters)[:-2] + ")"

    def store_force(self, dof, force):
        self.added_masses[dof] = force.real
        self.radiation_dampings[dof] = self.problem.omega * force.imag

    def records(self):
        return [dict(omega=self.omega, influenced_dof=dof, radiating_dof=self.radiating_dof,
                     added_mass=self.added_masses[dof], radiation_damping=self.radiation_dampings[dof])
                for dof in self.influenced_dofs]


def assemble_dataset(results):
    import pandas as pd
    import xarray as xr

    df = pd.DataFrame([record for result in results for record in result.records()])

    dataset = xr.Dataset()

    if 'added_mass' in df.columns:
        radiation_cases = df[df['added_mass'].notnull()].dropna(1)
        radiation_cases = radiation_cases.set_index(['omega', 'radiating_dof', 'influenced_dof'])
        radiation_cases = radiation_cases.to_xarray()
        dataset = xr.merge([dataset, radiation_cases])

    if 'diffraction_force' in df.columns:
        diffraction_cases = df[df['diffraction_force'].notnull()].dropna(1)
        diffraction_cases = diffraction_cases.set_index(['omega', 'angle', 'influenced_dof'])
        diffraction_cases = diffraction_cases.to_xarray()
        dataset = xr.merge([dataset, diffraction_cases])

    return dataset

