=======
Outputs
=======

Save the dataset as NetCDF file
-------------------------------

The xarray dataset produced by :func:`assemble_dataset <capytaine.results.assemble_dataset>` (or :meth:`fill_dataset <capytaine.Nemoh.Nemoh.fill_dataset>`) has a structure close to the NetCDF file format and can easily be saved to this format.
However, the netCDF standard does not handle complex numbers.
The complex-valued array can be saved as a bigger real-valued array with the helf of the :mod:`capytaine.io.xarray` module::

    from capytaine.io.xarray import separate_complex_values
    separate_complex_values(dataset).to_netcdf("path/to/dataset.nc")

The dataset can then be reloaded by::

    from capytaine.io.xarray import merge_complex_values
    dataset = merge_complex_values(xr.open_dataset("path/to/dataset.nc"))

See also the `documentation of xarray`_ for more details and options.

.. _`documentation of xarray`: http://xarray.pydata.org/en/stable/io.html

.. note:: If you get the error message 
   :code:`ValueError: unsupported dtype for netCDF4 variable: object`,
   please check `this Github issue <https://github.com/mancellin/capytaine/issues/2>`_.


