from itertools import product

import numpy as np
import pandas as pd
import xarray as xr

from capytaine.bem.problems_and_results import LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem


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


ATTRIBUTE_RATHER_THAN_COORD = True
# If True, the coordinates with a single value are replaced by attributes in the dataset.


def _squeeze_dimensions(data_array, dimensions=None):
    """Remove dimensions if they are of size 1."""
    if dimensions is None:
        dimensions = data_array.dims
    for dim in dimensions:
        if len(data_array[dim]) == 1:
            data_array = data_array.squeeze(dim, drop=ATTRIBUTE_RATHER_THAN_COORD)
    return data_array


def wavenumber_data_array(results):
    """Read the wavenumber in a list of :class:`LinearPotentialFlowResult`
    and store them into a :class:`xarray.DataArray`."""
    records = [dict(g=result.g, water_depth=result.depth, omega=result.omega, wavenumber=result.wavenumber)
               for result in results]
    optional_vars = ['g', 'water_depth']
    dimensions = ['omega']
    df = pd.DataFrame(records).drop_duplicates()
    df = df.set_index(optional_vars + dimensions)
    array = df.to_xarray()['wavenumber']
    array = _squeeze_dimensions(array, dimensions=optional_vars)
    return array


def assemble_dataset(results):
    """Transform a list of :class:`LinearPotentialFlowResult` to a :class:`xarray.Dataset`."""
    dataset = xr.Dataset()

    df = pd.DataFrame([record for result in results for record in result.records])
    if len(df) == 0:
        raise ValueError("No result passed to assemble_dataset.")

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
        diffraction_cases = diffraction_cases.set_index(optional_vars + ['convention'] + dimensions)
        diffraction_cases = diffraction_cases.to_xarray()
        diffraction_cases = _squeeze_dimensions(diffraction_cases, dimensions=optional_vars + ['convention'])
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
