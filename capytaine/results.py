#!/usr/bin/env python
# coding: utf-8
"""
Definition of the objects storing the results of a LinearPotentialFlowProblem.
"""

import logging

from attr import attrs, attrib, Factory
import numpy as np



LOG = logging.getLogger(__name__)


@attrs(slots=True)
class LinearPotentialFlowResult:
    problem = attrib()
    forces = attrib(default=Factory(dict), init=False, repr=False)

    sources = attrib(default=None, init=False, repr=False)
    potential = attrib(default=None, init=False, repr=False)

    def store_force(self, dof, force):
        self.forces[dof] = force

    def __getattr__(self, name):
        try:
            return getattr(self.problem, name)
        except AttributeError:
            raise AttributeError(f"{self.__class__} does not have a attribute named {name}.")


@attrs(slots=True)
class DiffractionResult(LinearPotentialFlowResult):
    pass


@attrs(slots=True)
class RadiationResult(LinearPotentialFlowResult):
    forces = attrib(default=None, init=False, repr=False)
    added_masses = attrib(default=Factory(dict), init=False, repr=False)
    radiation_dampings = attrib(default=Factory(dict), init=False, repr=False)

    def store_force(self, dof, force):
        self.added_masses[dof] = force.real
        self.radiation_dampings[dof] = self.problem.omega * force.imag


def assemble_radiation_results_matrices(results):
    """Combine the results of several radiation problems into an array of added masses and an array of radiation dampings.

    Parameters
    ----------
    results: list of RadiationResults
        The results container from which to extract the date.
        If objects other than RadiationResults are in the list, they are silently ignored.

    Returns
    -------
    added_masses: 3D xarray
        the added masses organised by frequencies and dofs
    radiation_dampings: 3D xarray
        the radiations_dampings organised by frequencies and dofs
    """

    import xarray as xr

    omegas = set()
    radiating_dofs = []
    influenced_dofs = []
    for result in results:
        if isinstance(result, RadiationResult):
            omegas.add(result.omega)
            if result.radiating_dof not in radiating_dofs:
                radiating_dofs.append(result.radiating_dof)
            for dof in result.influenced_dofs:
                if dof not in influenced_dofs:
                    influenced_dofs.append(dof)

    omegas = sorted(list(omegas))

    added_masses = xr.DataArray(np.empty((len(omegas), len(radiating_dofs), len(influenced_dofs))),
                                dims=('omega', 'radiating_dof', 'influenced_dof'),
                                coords={'omega': omegas, 'radiating_dof': radiating_dofs, 'influenced_dof': influenced_dofs})
    radiation_dampings = added_masses.copy(deep=True)

    for result in results:
        if isinstance(result, RadiationResult):
            for dof in result.influenced_dofs:
                added_masses.loc[dict(omega=result.omega,
                                      radiating_dof=result.radiating_dof,
                                      influenced_dof=dof)] = result.added_masses[dof]
                radiation_dampings.loc[dict(omega=result.omega,
                                            radiating_dof=result.radiating_dof,
                                            influenced_dof=dof)] = result.radiation_dampings[dof]

    return added_masses, radiation_dampings
