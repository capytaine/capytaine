#!/usr/bin/env python
# coding: utf-8
"""Definition of the objects storing the results of a LinearPotentialFlowProblem."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging

from attr import attrs, attrib, Factory

from capytaine.tools.Airy_wave import Froude_Krylov_force


LOG = logging.getLogger(__name__)


@attrs(slots=True)
class LinearPotentialFlowResult:
    problem = attrib()

    sources = attrib(default=None, init=False, repr=False)
    potential = attrib(default=None, init=False, repr=False)

    fs_elevation = attrib(default=Factory(dict), init=False, repr=False)

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

    @property
    def settings_dict(self):
        return dict(omega=self.omega, angle=self.angle,
                    water_depth=self.depth, body_name=self.body.name, g=self.g, rho=self.rho)

    @property
    def records(self):
        return [dict(self.settings_dict, influenced_dof=dof,
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

    @property
    def settings_dict(self):
            return dict(omega=self.omega, radiating_dof=self.radiating_dof,
                        water_depth=self.depth, body_name=self.body.name, g=self.g, rho=self.rho)

    @property
    def records(self):
        return [dict(self.settings_dict, influenced_dof=dof,
                     added_mass=self.added_masses[dof], radiation_damping=self.radiation_dampings[dof])
                for dof in self.influenced_dofs]


##############################################################################################
#                                     xarray assembling                                      #
##############################################################################################

ATTRIBUTE_RATHER_THAN_COORD = True
# If True, the coordinates with a single value are replaced by attributes in the dataset.


def _squeeze_dimensions(data_array, dimensions=None):
    if dimensions is None:
        dimensions = data_array.dims
    for dim in dimensions:
        if len(data_array[dim]) == 1:
            data_array = data_array.squeeze(dim, drop=ATTRIBUTE_RATHER_THAN_COORD)
    return data_array


def add_wavenumber_coord(dataset, results):
    import pandas as pd
    wavenumbers = [dict(g=result.g, water_depth=result.depth, omega=result.omega, wavenumber=result.wavenumber) for result in results]
    df = pd.DataFrame(wavenumbers).drop_duplicates()
    df = df.set_index(['g', 'water_depth', 'omega'])
    dataset.coords['wavenumber'] = (['g', 'water_depth', 'omega'], df.to_xarray()['wavenumber'].data)


def assemble_dataset(results):
    import pandas as pd
    import xarray as xr

    dataset = xr.Dataset()

    df = pd.DataFrame([record for result in results for record in result.records])

    optional_vars = ['g', 'rho', 'body_name', 'water_depth']

    # RADIATION RESULTS
    if 'added_mass' in df.columns:
        radiation_cases = df[df['added_mass'].notnull()].dropna(1)

        dimensions = ['omega', 'radiating_dof', 'influenced_dof']
        radiation_cases = radiation_cases.set_index(optional_vars + dimensions)
        radiation_cases = radiation_cases.to_xarray()
        radiation_cases = _squeeze_dimensions(radiation_cases, dimensions=optional_vars)
        dataset = xr.merge([dataset, radiation_cases])

    # DIFFRACTION RESULTS
    if 'diffraction_force' in df.columns:
        diffraction_cases = df[df['diffraction_force'].notnull()].dropna(1)

        dimensions = ['omega', 'angle', 'influenced_dof']
        diffraction_cases = diffraction_cases.set_index(optional_vars + dimensions)
        diffraction_cases = diffraction_cases.to_xarray()
        diffraction_cases = _squeeze_dimensions(diffraction_cases, dimensions=optional_vars)
        dataset = xr.merge([dataset, diffraction_cases])

    # BODIES PROPERTIES
    bodies = list({result.body for result in results})
    for body_property in ['mass', 'hydrostatic_stiffness']:
        bodies_properties = {body.name: body.__getattribute__(body_property) for body in bodies if hasattr(body, body_property)}
        if len(bodies_properties) > 0:
            bodies_properties = xr.concat(bodies_properties.values(), pd.Index(bodies_properties.keys(), name='body_name'))
            bodies_properties = _squeeze_dimensions(bodies_properties, dimensions=['body_name'])
            dataset = xr.merge([dataset, {body_property: bodies_properties}])

    # ATTRIBUTES
    if ATTRIBUTE_RATHER_THAN_COORD:
        for optional_var in optional_vars:
            optional_var_range = df[optional_var].unique()
            if len(optional_var_range) == 1:
                dataset.attrs[optional_var] = optional_var_range[0]

    return dataset

