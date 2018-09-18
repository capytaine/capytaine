
import numpy as np
import xarray as xr


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


