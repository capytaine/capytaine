import logging
from datetime import datetime
from itertools import product
from typing import Sequence

import numpy as np
import pandas as pd
import xarray as xr

from capytaine.bem.problems_and_results import \
    LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem, \
    LinearPotentialFlowResult
from capytaine.post_pro.kochin import compute_kochin


LOG = logging.getLogger(__name__)


def problems_from_dataset(dataset, bodies):
    """Generate a list of problems from the coordinates of a dataset.

    Parameters
    ----------
    dataset : xarray Dataset
        dataset containing the problems parameters: frequency, radiating_dof, water_depth, ...
    bodies : list of FloatingBody
        the bodies involved in the problems

    Returns
    -------
    list of LinearPotentialFlowProblem
    """
    assert len(list(set(body.name for body in bodies))) == len(bodies), \
        "All bodies should have different names."

    dataset = _unsqueeze_dimensions(dataset)

    omega_range = dataset['omega'].data if 'omega' in dataset else [LinearPotentialFlowProblem.default_parameters['omega']]
    angle_range = dataset['angle'].data if 'angle' in dataset else None
    radiating_dofs = dataset['radiating_dof'].data if 'radiating_dof' in dataset else None
    water_depth_range = dataset['water_depth'].data if 'water_depth' in dataset else [-LinearPotentialFlowProblem.default_parameters['sea_bottom']]
    rho_range = dataset['rho'].data if 'rho' in dataset else [LinearPotentialFlowProblem.default_parameters['rho']]

    if 'body_name' in dataset:
        assert set(dataset['body_name'].data) <= {body.name for body in bodies}
        body_range = {body.name: body for body in bodies if body.name in dataset['body_name'].data}
    else:
        body_range = {body.name: body for body in bodies}

    problems = []
    if angle_range is not None:
        for omega, angle, water_depth, body_name, rho \
                in product(omega_range, angle_range, water_depth_range, body_range, rho_range):
            problems.append(
                DiffractionProblem(body=body_range[body_name], omega=omega,
                                   angle=angle, sea_bottom=-water_depth, rho=rho)
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
    """Remove dimensions if they are of size 1."""
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


def _dataset_from_dataframe(df: pd.DataFrame, variables, dimensions, optional_dims):
    for variable_name in variables:
        df = df[df[variable_name].notnull()].dropna(1)  # Keep only records with non null values of all the variables
    df = df.drop_duplicates()
    df = df.set_index(optional_dims + dimensions)

    da = df.to_xarray()[variables]
    da = _squeeze_dimensions(da, dimensions=optional_dims)
    return da


def wavenumber_data_array(results: Sequence[LinearPotentialFlowResult]):
    """Read the wavenumber in a list of :class:`LinearPotentialFlowResult`
    and store them into a :class:`xarray.DataArray`."""
    records = pd.DataFrame(
        [dict(g=result.g, water_depth=result.depth, omega=result.omega, wavenumber=result.wavenumber)
         for result in results]
    )
    ds = _dataset_from_dataframe(records, variables=['wavenumber'], dimensions=['omega'], optional_dims=['g', 'water_depth'])
    return ds['wavenumber']


def hydrostatics_dataset(bodies):
    dataset = xr.Dataset()
    for body_property in ['mass', 'hydrostatic_stiffness']:
        bodies_properties = {body.name: body.__getattribute__(body_property) for body in bodies if hasattr(body, body_property)}
        if len(bodies_properties) > 0:
            bodies_properties = xr.concat(bodies_properties.values(), pd.Index(bodies_properties.keys(), name='body_name'))
            bodies_properties = _squeeze_dimensions(bodies_properties, dimensions=['body_name'])
            dataset = xr.merge([dataset, {body_property: bodies_properties}])
    return dataset


def kochin_data_array(results, theta_range, **kwargs):
    records = pd.DataFrame([dict(result.settings_dict, theta=theta, kochin=kochin)
                            for result in results
                            for theta, kochin in zip(theta_range, compute_kochin(result, theta_range, **kwargs))])

    ds = _dataset_from_dataframe(records, ['kochin'],
                                 dimensions=['omega', 'radiating_dof', 'theta'],
                                 optional_dims=['g', 'rho', 'body_name', 'water_depth'])
    return ds['kochin']


def assemble_dataset(results: Sequence[LinearPotentialFlowResult],
                     wavenumber=False, wavelength=False, mesh=False, hydrostatics=True,
                     attrs=None):
    """Transform a list of :class:`LinearPotentialFlowResult` to a :class:`xarray.Dataset`."""
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
            dimensions=['omega', 'angle', 'influenced_dof'],
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
        dataset.coords['nb_faces'] = ('body_name', [nb_faces])

    # HYDROSTATICS
    if hydrostatics:
        bodies = list({result.body for result in results})
        dataset = xr.merge([dataset, hydrostatics_dataset(bodies)])

    dataset.attrs.update(attrs)
    return dataset


def separate_complex_values(ds: xr.Dataset) -> xr.Dataset:
    """Return a new Dataset where complex-valued arrays of shape (...) have been replaced by real-valued arrays of shape (2, ...).
    Invert of :func:`merge_complex_values`."""
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
    """Return a new Dataset where real-valued arrays of shape (2, ...) have been replaced by complex-valued arrays of shape (...).
    Invert of :func:`separate_complex_values`."""
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


def save_in_netcdf(ds: xr.Dataset, filepath):
    ds.to_netcdf(filepath,
                 encoding={'radiating_dof': {'dtype': 'U'},
                           'influenced_dof': {'dtype': 'U'}})
