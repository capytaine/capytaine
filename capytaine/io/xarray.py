"""Tools to use xarray Datasets as inputs and outputs.

.. todo:: This module could be tidied up a bit and some methods merged or
          uniformized.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from datetime import datetime
from itertools import product
from collections import Counter
from typing import Sequence, List, Union

import numpy as np
import pandas as pd
import xarray as xr

from capytaine import __version__
from capytaine.bodies.bodies import FloatingBody
from capytaine.bem.problems_and_results import (
    LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem,
    LinearPotentialFlowResult, _default_parameters)
from capytaine.post_pro.kochin import compute_kochin
from capytaine.io.bemio import dataframe_from_bemio


LOG = logging.getLogger(__name__)


#########################
#  Reading test matrix  #
#########################

def problems_from_dataset(dataset: xr.Dataset,
                          bodies: Union[FloatingBody, Sequence[FloatingBody]],
                          ) -> List[LinearPotentialFlowProblem]:
    """Generate a list of problems from a test matrix.

    Parameters
    ----------
    dataset : xarray Dataset
        Test matrix containing the problems parameters.
    bodies : FloatingBody or list of FloatingBody
        The bodies on which the computations of the test matrix will be applied.
        They should all have different names.

    Returns
    -------
    list of LinearPotentialFlowProblem

    Raises
    ------
    ValueError
        if required fields are missing in the dataset
    """
    if isinstance(bodies, FloatingBody):
        bodies = [bodies]

    # Should be done before looking for `frequency_keys`, otherwise
    # frequencies provided as a scalar dimension will be skipped.
    dataset = _unsqueeze_dimensions(dataset)

    # SANITY CHECKS
    assert len(list(set(body.name for body in bodies))) == len(bodies), \
        "All bodies should have different names."

    # Warn user in case of key with unrecognized name (e.g. misspells)
    keys_in_dataset = set(dataset.dims)
    accepted_keys = {'wave_direction', 'radiating_dof', 'influenced_dof',
                     'body_name', 'omega', 'period', 'wavelength', 'wavenumber',
                     'forward_speed', 'water_depth', 'rho', 'g', 'theta'}
    unrecognized_keys = keys_in_dataset.difference(accepted_keys)
    if len(unrecognized_keys) > 0:
        LOG.warning(f"Unrecognized key(s) in dataset: {unrecognized_keys}")

    if ("radiating_dof" not in keys_in_dataset) and ("wave_direction" not in keys_in_dataset):
        raise ValueError("Neither 'radiating_dof' nor 'wave_direction' has been provided in the dataset. "
                "No linear potential flow problem can be inferred.")

    frequency_keys = keys_in_dataset & {'omega', 'period', 'wavelength', 'wavenumber'}
    if len(frequency_keys) > 1:
            raise ValueError("Setting problems requires at most one of the following: omega (angular frequency) OR period OR wavenumber OR wavelength.\n"
                             "Received {}".format(frequency_keys))
    # END SANITY CHECKS

    if len(frequency_keys) == 0:
        freq_type = "omega"
        freq_range = [_default_parameters['omega']]
    else:  # len(frequency_keys) == 1
        freq_type = list(frequency_keys)[0]  # Get the only item
        freq_range = dataset[freq_type].data

    water_depth_range = dataset['water_depth'].data if 'water_depth' in dataset else [_default_parameters['water_depth']]
    rho_range = dataset['rho'].data if 'rho' in dataset else [_default_parameters['rho']]
    g_range = dataset['g'].data if 'g' in dataset else [_default_parameters['g']]
    forward_speed_range = dataset['forward_speed'] if 'forward_speed' in dataset else [_default_parameters['forward_speed']]

    wave_direction_range = dataset['wave_direction'].data if 'wave_direction' in dataset else None
    radiating_dofs = dataset['radiating_dof'].data.astype(object) if 'radiating_dof' in dataset else None
    # astype(object) is meant to convert Numpy internal string type numpy.str_ to Python general string type.

    if 'body_name' in dataset:
        assert set(dataset['body_name'].data) <= {body.name for body in bodies}, \
            "Some body named in the dataset was not given as argument to `problems_from_dataset`."
        body_range = {body.name: body for body in bodies if body.name in dataset['body_name'].data}
        # Only the bodies listed in the dataset have been kept
    else:
        body_range = {body.name: body for body in bodies}

    problems = []
    if wave_direction_range is not None:
        for freq, wave_direction, water_depth, body_name, forward_speed, rho, g \
                in product(freq_range, wave_direction_range, water_depth_range, body_range, forward_speed_range, rho_range, g_range):
            if freq not in {0.0, np.inf}:
                problems.append(
                    DiffractionProblem(body=body_range[body_name], **{freq_type: freq},
                                       wave_direction=wave_direction, water_depth=water_depth,
                                       forward_speed=forward_speed, rho=rho, g=g)
                )
            elif freq in {0.0, np.inf} and radiating_dofs is not None:
                # Diffraction problems are not defined for 0 and infinite frequency.
                # But we don't want the whole batch to fail, as these frequencies are there for the radiation problems.
                # The excitation force will be NaN for these frequencies in the resulting dataset.
                pass
            else:
                raise ValueError("Zero and infinite frequencies are not defined when solving only diffraction problems.")

    if radiating_dofs is not None:
        for freq, radiating_dof, water_depth, body_name, forward_speed, rho, g \
                in product(freq_range, radiating_dofs, water_depth_range, body_range, forward_speed_range, rho_range, g_range):
            if forward_speed == 0.0:
                problems.append(
                    RadiationProblem(body=body_range[body_name], **{freq_type: freq},
                                     radiating_dof=radiating_dof, water_depth=water_depth,
                                     forward_speed=forward_speed, rho=rho, g=g)
                )
            else:
                if wave_direction_range is None:
                    LOG.warning("Dataset contains non-zero forward speed (forward_speed=%.2f) but no wave_direction has been provided. Wave direction of 0 rad (x-axis) has been assumed.", forward_speed)
                    wave_direction_range = [0.0]
                for wave_direction in wave_direction_range:
                    problems.append(
                        RadiationProblem(body=body_range[body_name], **{freq_type: freq},
                                         radiating_dof=radiating_dof, water_depth=water_depth,
                                         forward_speed=forward_speed, wave_direction=wave_direction,
                                         rho=rho, g=g)
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
        df = df[df[variable_name].notnull()].dropna(axis='columns')  # Keep only records with non null values of all the variables
    df = df.drop_duplicates(optional_dims + dimensions)
    df = df.set_index(optional_dims + dimensions)

    da = df.to_xarray()[variables]
    da = _squeeze_dimensions(da, dimensions=optional_dims)
    return da


def hydrostatics_dataset(bodies: Sequence[FloatingBody]) -> xr.Dataset:
    """Create a dataset by looking for 'inertia_matrix' and 'hydrostatic_stiffness'
    for each of the bodies in the list passed as argument.
    """
    dataset = xr.Dataset()
    for body_property in ['inertia_matrix', 'hydrostatic_stiffness']:
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
    records = pd.DataFrame([
        dict(**result.problem._asdict(), theta=theta, kochin=kochin, kind=result.__class__.__name__)
        for result in results
        for theta, kochin in zip(theta_range.data,
                                 compute_kochin(result, theta_range, **kwargs))
    ])

    kochin_data = xr.Dataset()

    if "RadiationResult" in set(records['kind']):
        radiation = _dataset_from_dataframe(
            records[records['kind'] == "RadiationResult"],
            variables=['kochin'],
            dimensions=['omega', 'radiating_dof', 'theta'],
            optional_dims=['g', 'rho', 'body_name', 'water_depth', 'forward_speed', 'wave_direction']
        )
        kochin_data['kochin'] = radiation['kochin']

    if "DiffractionResult" in set(records['kind']):
        diffraction = _dataset_from_dataframe(
            records[records['kind'] == "DiffractionResult"],
            ['kochin'],
            dimensions=['omega', 'wave_direction', 'theta'],
            optional_dims=['g', 'rho', 'body_name', 'water_depth', 'forward_speed']
        )
        kochin_data['kochin_diffraction'] = diffraction['kochin']

    return kochin_data


def collect_records(results):
    records_list = []
    warned_once_about_no_free_surface = False
    for result in results:
        if result.free_surface == np.inf:
            if not warned_once_about_no_free_surface:
                LOG.warning("Datasets currently only support cases with a free surface (free_surface=0.0).\n"
                            "Cases without a free surface (free_surface=inf) are ignored.\n"
                            "See also https://github.com/mancellin/capytaine/issues/88")
                warned_once_about_no_free_surface = True
            else:
                pass
        else:
            for record in result.records:
                records_list.append(record)
    return records_list

def assemble_dataset(results,
                     omega=True, wavenumber=True, wavelength=True, period=True,
                     mesh=False, hydrostatics=True, attrs=None) -> xr.Dataset:
    """Transform a list of :class:`LinearPotentialFlowResult` into a :class:`xarray.Dataset`.

    .. todo:: The :code:`mesh` option to store information on the mesh could be improved.
              It could store the full mesh in the dataset to ensure the reproducibility of
              the results.

    Parameters
    ----------
    results: list of LinearPotentialFlowResult
        The results that will be read.
    omega: bool, optional
        If True, the coordinate 'omega' will be added to the output dataset.
    wavenumber: bool, optional
        If True, the coordinate 'wavenumber' will be added to the output dataset.
    wavelength: bool, optional
        If True, the coordinate 'wavelength' will be added to the output dataset.
    period: bool, optional
        If True, the coordinate 'period' will be added to the output dataset.
    mesh: bool, optional
        If True, store some infos on the mesh in the output dataset.
    hydrostatics: bool, optional
        If True, store the hydrostatic data in the output dataset if they exist.
    attrs: dict, optional
        Attributes that should be added to the output dataset.
    """
    dataset = xr.Dataset()

    error_msg = 'The first argument of `assemble_dataset` must be either a list of LinearPotentialFlowResult or a bemio.io object'
    if hasattr(results, '__iter__'):
        try:
            if 'capytaine' in results[0].__module__:
                bemio_import = False
            else:
                raise TypeError(error_msg)
        except:
            raise TypeError(error_msg)

    else:
        try:
            if 'bemio.io' in results.__module__:
                bemio_import = True
            else:
                raise TypeError(error_msg)
        except:
            raise TypeError(error_msg)

    if bemio_import:
        records = dataframe_from_bemio(results, wavenumber, wavelength) # TODO add hydrostatics
        all_dofs_in_order = {'Surge': None, 'Sway': None, 'Heave': None, 'Roll': None, 'Pitch': None, 'Yaw': None}
        main_freq_type = "omega"

    else:
        records = pd.DataFrame(collect_records(results))
        all_dofs_in_order = {k: None for r in results for k in r.body.dofs.keys()}
        main_freq_type = Counter((res.provided_freq_type for res in results)).most_common(1)[0][0]

    if attrs is None:
        attrs = {}
    attrs['creation_of_dataset'] = datetime.now().isoformat()

    if len(records) == 0:
        raise ValueError("No result passed to assemble_dataset.")

    inf_dof_cat = pd.CategoricalDtype(categories=all_dofs_in_order.keys())
    records["influenced_dof"] = records["influenced_dof"].astype(inf_dof_cat)
    rad_dof_cat = pd.CategoricalDtype(categories=all_dofs_in_order.keys())
    if 'added_mass' in records.columns:
        records["radiating_dof"] = records["radiating_dof"].astype(rad_dof_cat)

    optional_dims = ['g', 'rho', 'body_name', 'water_depth', 'forward_speed']

    # RADIATION RESULTS
    if 'added_mass' in records.columns:
        radiation_cases = _dataset_from_dataframe(
            records,
            variables=['added_mass', 'radiation_damping'],
            dimensions=[main_freq_type, 'radiating_dof', 'influenced_dof'],
            optional_dims=optional_dims + ['wave_direction'])
        radiation_cases.added_mass.attrs['long_name'] = 'Added mass'
        radiation_cases.radiation_damping.attrs['long_name'] = 'Radiation damping'
        radiation_cases.radiating_dof.attrs['long_name'] = 'Radiating DOF'
        radiation_cases.influenced_dof.attrs['long_name'] = 'Influenced DOF'
        dataset = xr.merge([dataset, radiation_cases])

    # DIFFRACTION RESULTS
    if 'diffraction_force' in records.columns:
        diffraction_cases = _dataset_from_dataframe(
            records,
            variables=['diffraction_force', 'Froude_Krylov_force'],
            dimensions=[main_freq_type, 'wave_direction', 'influenced_dof'],
            optional_dims=optional_dims)
        diffraction_cases.diffraction_force.attrs['long_name'] = 'Diffraction force'
        diffraction_cases.Froude_Krylov_force.attrs['long_name'] = 'Froude Krylov force'
        diffraction_cases.influenced_dof.attrs['long_name'] = 'Influenced DOF'
        diffraction_cases.wave_direction.attrs['long_name'] = 'Wave direction'
        diffraction_cases.wave_direction.attrs['units'] = 'rad'
        dataset = xr.merge([dataset, diffraction_cases])
        dataset['excitation_force'] = dataset['Froude_Krylov_force'] + dataset['diffraction_force']

    # OTHER FREQUENCIES TYPES
    if omega and main_freq_type != "omega":
        omega_ds = _dataset_from_dataframe(
                records,
                variables=['omega'],
                dimensions=[main_freq_type],
                optional_dims=['g', 'water_depth'] if main_freq_type in {'wavelength', 'wavenumber'} else []
                )
        dataset.coords['omega'] = omega_ds['omega']
        dataset.omega.attrs['long_name'] = 'Angular frequency'
        dataset.omega.attrs['units'] = 'rad/s'

    if period and main_freq_type != "period":
        period_ds = _dataset_from_dataframe(
                records,
                variables=['period'],
                dimensions=[main_freq_type],
                optional_dims=['g', 'water_depth'] if main_freq_type in {'wavelength', 'wavenumber'} else []
                )
        dataset.coords['period'] = period_ds['period']
        dataset.period.attrs['long_name'] = 'Period'
        dataset.period.attrs['units'] = 's'

    if wavenumber and main_freq_type != "wavenumber":
        wavenumber_ds = _dataset_from_dataframe(
                records,
                variables=['wavenumber'],
                dimensions=[main_freq_type],
                optional_dims=['g', 'water_depth'] if main_freq_type in {'period', 'omega'} else []
                )
        dataset.coords['wavenumber'] = wavenumber_ds['wavenumber']
        dataset.wavenumber.attrs['long_name'] = 'Angular wavenumber'
        dataset.wavenumber.attrs['units'] = 'rad/m'

    if wavelength and main_freq_type != "wavelength":
        wavelength_ds = _dataset_from_dataframe(
                records,
                variables=['wavelength'],
                dimensions=[main_freq_type],
                optional_dims=['g', 'water_depth'] if main_freq_type in {'period', 'omega'} else []
                )
        dataset.coords['wavelength'] = wavelength_ds['wavelength']
        dataset.wavelength.attrs['long_name'] = 'Wave length'
        dataset.wavelength.attrs['units'] = 'm'

    if not all(records["forward_speed"] == 0.0):
        omegae_ds = _dataset_from_dataframe(
                records,
                variables=['encounter_omega'],
                dimensions=['forward_speed', 'wave_direction', main_freq_type],
                optional_dims=['g', 'water_depth'],
                )
        dataset.coords['encounter_omega'] = omegae_ds['encounter_omega']
        dataset.encounter_omega.attrs['long_name'] = 'Encounter angular frequency'
        dataset.encounter_omega.attrs['units'] = 'rad/s'

        encounter_wave_direction_ds = _dataset_from_dataframe(
                records,
                variables=['encounter_wave_direction'],
                dimensions=['forward_speed', 'wave_direction', main_freq_type],
                optional_dims=[],
                )
        dataset.coords['encounter_wave_direction'] = encounter_wave_direction_ds['encounter_wave_direction']
        dataset.encounter_wave_direction.attrs['long_name'] = 'Encounter wave direction'
        dataset.encounter_wave_direction.attrs['units'] = 'rad'

    if mesh:
        if bemio_import:
            LOG.warning('Bemio data does not include mesh data. mesh=True is ignored.')
        else:
            # TODO: Store full mesh...
            bodies = list({result.body for result in results})  # Filter out duplicate bodies in the list of results
            nb_faces = {body.name: body.mesh.nb_faces for body in bodies}

            def name_or_str(c):
                return c.name if hasattr(c, 'name') else str(c)
            quad_methods = {body.name: name_or_str(body.mesh.quadrature_method) for body in bodies}

            if len(nb_faces) > 1:
                dataset.coords['nb_faces'] = ('body_name', [nb_faces[name] for name in dataset.coords['body_name'].data])
                dataset.coords['quadrature_method'] = ('body_name', [quad_methods[name] for name in dataset.coords['body_name'].data])
            else:
                def the_only(d):
                    """Return the only element of a 1-element dictionary"""
                    return next(iter(d.values()))
                dataset.coords['nb_faces'] = the_only(nb_faces)
                dataset.coords['quadrature_method'] = the_only(quad_methods)

    # HYDROSTATICS
    if hydrostatics:
        if bemio_import:
            LOG.warning('Bemio data import being used, hydrostatics=True is ignored.')
        else:
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
        if ds[variable].dtype == complex:
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
        ds = ds.drop_vars('complex')
    return ds
