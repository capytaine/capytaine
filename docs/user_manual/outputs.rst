=======
Outputs
=======

Save the dataset as NetCDF file
-------------------------------

The xarray dataset produced by :code:`assembre_results` has a structure close to the NetCDF file format.
It can be saved in this format using for instance the following syntax::

    dataset.to_netcdf("path/to/dataset.nc", engine="h5netcdf")

The dataset can be reloaded by::

    dataset = xr.open_dataset("path/to/dataset.nc", engine="h5netcdf")

See also the `documentation of xarray`_ for more details and options.
Note that the code above requires the package :code:`h5netcdf`, which can be installed with::

    $ pip install h5netcdf

.. _`documentation of xarray`: http://xarray.pydata.org/en/stable/io.html

.. warning:: The NetCDF standard does not handle complex numbers.
    The above code is able to save complex data but it might not be compatible with other software.
    In the future, a more robust implementation should be added to Capytaine.

