#!/usr/bin/env python
# coding: utf-8
"""Tools to use xarray Datasets as inputs and outputs.

.. todo:: This module could be tidied up a bit and some methods merged or
          uniformized.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from datetime import datetime
from itertools import product
from typing import Sequence, List, Union

import numpy as np
import pandas as pd
import xarray as xr

from capytaine import __version__
from capytaine.bodies.bodies import FloatingBody
from capytaine.bem.problems_and_results import (
    LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem,
    LinearPotentialFlowResult)
from capytaine.post_pro.kochin import compute_kochin


LOG = logging.getLogger(__name__)


#########################
#  Reading test matrix  #
#########################

def problems_from_dataset(dataset: xr.Dataset,
                          bodies: Sequence[FloatingBody],
                          ) -> List[LinearPotentialFlowProblem]:
    """Generate a list of problems from a test matrix.

    Parameters
    ----------
    dataset : xarray Dataset
        Test matrix containing the problems parameters.
    bodies : list of FloatingBody
        The bodies on which the computations of the test matrix will be applied.
        They should all have different names.

    Returns
    -------
    list of LinearPotentialFlowProblem
    """
    assert len(list(set(body.name for body in bodies))) == len(bodies), \
        "All bodies should have different names."

    dataset = _unsqueeze_dimensions(dataset)

    omega_range = dataset['omega'].data if 'omega' in dataset else [LinearPotentialFlowProblem.default_parameters['omega']]
    water_depth_range = dataset['water_depth'].data if 'water_depth' in dataset else [-LinearPotentialFlowProblem.default_parameters['sea_bottom']]
    rho_range = dataset['rho'].data if 'rho' in dataset else [LinearPotentialFlowProblem.default_parameters['rho']]

    wave_direction_range = dataset['wave_direction'].data if 'wave_direction' in dataset else None
    radiating_dofs = dataset['radiating_dof'].data if 'radiating_dof' in dataset else None

    if 'body_name' in dataset:
        assert set(dataset['body_name'].data) <= {body.name for body in bodies}, \
            "Some body named in the dataset was not given as argument to `problems_from_dataset`."
        body_range = {body.name: body for body in bodies if body.name in dataset['body_name'].data}
        # Only the bodies listed in the dataset have been kept
    else:
        body_range = {body.name: body for body in bodies}

    problems = []
    if wave_direction_range is not None:
        for omega, wave_direction, water_depth, body_name, rho \
                in product(omega_range, wave_direction_range, water_depth_range, body_range, rho_range):
            problems.append(
                DiffractionProblem(body=body_range[body_name], omega=omega,
                                   wave_direction=wave_direction, sea_bottom=-water_depth, rho=rho)
            )

    if radiating_dofs is not None:
        for omega, radiating_dof, water_depth, body_name, rho \
                in product(omega_range, radiating_dofs, water_depth_range, body_range, rho_range):
            problems.append(
                RadiationProblem(body=body_range[body_name], omega=omega,
                                 radiating_dof=radiating_dof, sea_bottom=-water_depth, rho=rho)
            )

    return sorted(problems)


def _squeeze_dimensions(data_array, dimensions=None):
    """Remove dimensions if they are of size 1. The coordinates become scalar coordinates."""
    if dimensions is None:
        dimensions = data_array.dims
    for dim in dimensions:
        if len(data_array[dim]) == 1:
            data_array = data_array.squeeze(dim, drop=False)
    return data_array


def _unsqueeze_dimensions(data_array, dimensions=None):
    """Add scalar coordinates as dimensions of size 1."""
    if dimensions is None:
        dimensions = list(data_array.coords.keys())
    for dim in dimensions:
        if len(data_array.coords[dim].values.shape) == 0:
            data_array = xr.concat([data_array], dim=dim)
    return data_array


######################
#  Dataset creation  #
######################

def _dataset_from_dataframe(df: pd.DataFrame,
                            variables: Union[str, Sequence[str]],
                            dimensions: Sequence[str],
                            optional_dims: Sequence[str],
                            ) -> Union[xr.DataArray, xr.Dataset]:
    """Transform a pandas.Dataframe into a xarray.Dataset.

    Parameters
    ----------
    df: pandas.DataFrame
        the input dataframe
    variables: string or sequence of strings
        the variables that will be stored in the output dataset.
        If a single name is provided, a DataArray of this variable will be provided instead.
    dimensions: sequence of strings
        Names of dimensions the variables depends on.
        They will always appear as dimension in the output dataset.
    optional_dims: sequence of strings
        Names of dimensions the variables depends on.
        They will appears as dimension in the output dataset only if they have
        more than one different values.
    """
    for variable_name in variables:
        df = df[df[variable_name].notnull()].dropna(1)  # Keep only records with non null values of all the variables
    df = df.drop_duplicates()
    df = df.set_index(optional_dims + dimensions)

    da = df.to_xarray()[variables]
    da = _squeeze_dimensions(da, dimensions=optional_dims)
    return da


def wavenumber_data_array(results: Sequence[LinearPotentialFlowResult]) -> xr.DataArray:
    """Read the wavenumbers in a list of :class:`LinearPotentialFlowResult`
    and store them into a :class:`xarray.DataArray`.
    """
    records = pd.DataFrame(
        [dict(g=result.g, water_depth=result.depth, omega=result.omega, wavenumber=result.wavenumber)
         for result in results]
    )
    ds = _dataset_from_dataframe(records, variables=['wavenumber'], dimensions=['omega'], optional_dims=['g', 'water_depth'])
    return ds['wavenumber']


def hydrostatics_dataset(bodies: Sequence[FloatingBody]) -> xr.Dataset:
    """Create a dataset by looking for 'mass' and 'hydrostatic_stiffness'
    for each of the bodies in the list passed as argument.
    """
    dataset = xr.Dataset()
    for body_property in ['mass', 'hydrostatic_stiffness']:
        bodies_properties = {body.name: body.__getattribute__(body_property) for body in bodies if hasattr(body, body_property)}
        if len(bodies_properties) > 0:
            bodies_properties = xr.concat(bodies_properties.values(), pd.Index(bodies_properties.keys(), name='body_name'))
            bodies_properties = _squeeze_dimensions(bodies_properties, dimensions=['body_name'])
            dataset = xr.merge([dataset, {body_property: bodies_properties}])
    return dataset


def kochin_data_array(results: Sequence[LinearPotentialFlowResult],
                      theta_range: Sequence[float],
                      **kwargs,
                      ) -> xr.Dataset:
    """Compute the Kochin function for a list of results and fills a dataset.

    .. seealso::
        :meth:`~capytaine.post_pro.kochin.compute_kochin`
            The present function is just a wrapper around :code:`compute_kochin`.
    """
    records = pd.DataFrame([dict(result.settings_dict, theta=theta, kochin=kochin)
                            for result in results
                            for theta, kochin in zip(theta_range.data, compute_kochin(result, theta_range, **kwargs))])

    ds = _dataset_from_dataframe(records, ['kochin'],
                                 dimensions=['omega', 'radiating_dof', 'theta'],
                                 optional_dims=['g', 'rho', 'body_name', 'water_depth'])
    return ds['kochin']


def assemble_dataset(results: Sequence[LinearPotentialFlowResult],
                     wavenumber=False, wavelength=False, mesh=False, hydrostatics=True,
                     attrs=None) -> xr.Dataset:
    """Transform a list of :class:`LinearPotentialFlowResult` into a :class:`xarray.Dataset`.

    .. todo:: The :code:`mesh` option to store informations on the mesh could be improved.
              It could store the full mesh in the dataset to ensure the reproducibility of
              the results.

    Parameters
    ----------
    results: list of LinearPotentialFlowResult
        The results that will be read.
    wavenumber: bool, optional
        If True, the coordinate 'wavenumber' will be added to the ouput dataset.
    wavelength: bool, optional
        If True, the coordinate 'wavelength' will be added to the ouput dataset.
    mesh: bool, optional
        If True, store some infos on the mesh in the output dataset.
    hydrostatics: bool, optional
        If True, store the hydrostatic data in the output dataset if they exist.
    attrs: dict, optional
        Attributes that should be added to the output dataset.
    """
    dataset = xr.Dataset()

    if attrs is None:
        attrs = {}
    attrs['creation_of_dataset'] = datetime.now().isoformat()

    records = pd.DataFrame([record for result in results for record in result.records])
    if len(records) == 0:
        raise ValueError("No result passed to assemble_dataset.")

    optional_dims = ['g', 'rho', 'body_name', 'water_depth']

    # RADIATION RESULTS
    if 'added_mass' in records.columns:
        radiation_cases = _dataset_from_dataframe(
            records,
            variables=['added_mass', 'radiation_damping'],
            dimensions=['omega', 'radiating_dof', 'influenced_dof'],
            optional_dims=optional_dims)
        dataset = xr.merge([dataset, radiation_cases])

    # DIFFRACTION RESULTS
    if 'diffraction_force' in records.columns:
        conventions = set(records['convention'].dropna())
        if len(conventions) > 1:
            LOG.warning("Assembling a dataset mixing several conventions.")
        else:
            attrs['incoming_waves_convention'] = conventions.pop()

        diffraction_cases = _dataset_from_dataframe(
            records,
            variables=['diffraction_force', 'Froude_Krylov_force'],
            dimensions=['omega', 'wave_direction', 'influenced_dof'],
            optional_dims=optional_dims)
        dataset = xr.merge([dataset, diffraction_cases])

    # WAVENUMBER
    if wavenumber:
        dataset.coords['wavenumber'] = wavenumber_data_array(results)

    if wavelength:
        dataset.coords['wavelength'] = 2*np.pi/wavenumber_data_array(results)

    if mesh:
        # TODO: Store full mesh...
        bodies = list({result.body for result in results})
        nb_faces = {body.name: body.mesh.nb_faces for body in bodies}
        if len(nb_faces) > 1:
            dataset.coords['nb_faces'] = ('body_name', nb_faces)
        else:
            dataset.coords['nb_faces'] = list(nb_faces.items())[0][1]

    # HYDROSTATICS
    if hydrostatics:
        bodies = list({result.body for result in results})
        dataset = xr.merge([dataset, hydrostatics_dataset(bodies)])

    dataset.attrs.update(attrs)
    dataset.attrs['capytaine_version'] = __version__
    return dataset


################################
#  Handling of complex values  #
################################

def separate_complex_values(ds: xr.Dataset) -> xr.Dataset:
    """Return a new Dataset where complex-valued arrays of shape (...)
    have been replaced by real-valued arrays of shape (2, ...).

    .. seealso::
        :func:`merge_complex_values`
            The invert operation
    """
    ds = ds.copy()
    for variable in ds.data_vars:
        if ds[variable].dtype == np.complex:
            da = ds[variable]
            new_da = xr.DataArray(np.asarray((np.real(da).data, np.imag(da).data)),
                                  dims=('complex',) + da.dims)
            ds[variable] = new_da
            ds.coords['complex'] = ['re', 'im']
    return ds


def merge_complex_values(ds: xr.Dataset) -> xr.Dataset:
    """Return a new Dataset where real-valued arrays of shape (2, ...)
    have been replaced by complex-valued arrays of shape (...).

    .. seealso::
        :func:`separate_complex_values`
            The invert operation
    """
    if 'complex' in ds.coords:
        ds = ds.copy()
        for variable in ds.data_vars:
            if 'complex' in ds[variable].coords:
                da = ds[variable]
                new_dims = [d for d in da.dims if d != 'complex']
                new_da = xr.DataArray(da.sel(complex='re').data + 1j*da.sel(complex='im').data, dims=new_dims)
                ds[variable] = new_da
        ds = ds.drop('complex')
    return ds
