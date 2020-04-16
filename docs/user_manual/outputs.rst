=======
Outputs
=======

Building a dataset
------------------

If you have a list of :code:`LinearPotentialFlowResult`, you can assemble
them in a xarray dataset for a more convenient post-processing. Use the
following syntax::

   from capytaine import assemble_dataset
   dataset = assemble_dataset(list_of_results)

If you gave a test matrix to the :code:`BEMSolver.fill_dataset` method, the
output will directly be an xarray dataset.

Both :code:`assemble_dataset` and :code:`fill_dataset` accept some optional keyword
arguments to store more informations in the dataset:

- :code:`wavenumber` (default: :code:`False`): add the wavenumber of the
  incoming waves in the dataset.
- :code:`wavelength` (default: :code:`False`): add the wavelength of the
  incoming waves in the dataset.
- :code:`mesh` (default: :code:`False`): add some informations about the mesh in
  the dataset (number of faces, quadrature method).
- :code:`hydrostatics` (default: :code:`True`): if hydrostatics data are
  available in the :code:`FloatingBody`, they are added to the dataset. 

.. note:: The code does its best to keep the degrees of freedom in the same
          order as they have been provided by the user, but there is no
          guarantee that this will always be the case.
          Nonetheless, in all cases the labels will always match the data.
          So selecting a value by the name of the dof will always return the right one::

              data.sel(radiating_dof="Heave", influenced_dof="Heave")

          You can also manually reorder the dofs with the following syntax::

              sorted_dofs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]
              print(data.sel(radiating_dof=sorted_dofs, influenced_dof=sorted_dofs))


Saving the dataset as NetCDF file
---------------------------------

The xarray dataset produced by :func:`assemble_dataset <capytaine.results.assemble_dataset>` (or :meth:`fill_dataset <capytaine.bem.solver.BEMSolver.fill_dataset>`) has a structure close to the NetCDF file format and can easily be saved to this format::

	dataset.to_netcdf("path/to/dataset.nc")

See the `documentation of xarray <http://xarray.pydata.org/en/stable/io.html>`_ for details and options.

There are however a couple of issues you should be aware of:


Complex numbers
~~~~~~~~~~~~~~~

The netCDF standard does not handle complex numbers.
As a workaround, the complex-valued array can be saved as a bigger real-valued array with the help of the :mod:`capytaine.io.xarray` module::

    from capytaine.io.xarray import separate_complex_values
    separate_complex_values(dataset).to_netcdf("path/to/dataset.nc")

The dataset can then be reloaded by::

    from capytaine.io.xarray import merge_complex_values
    dataset = merge_complex_values(xr.open_dataset("path/to/dataset.nc"))


String format
~~~~~~~~~~~~~

There is an issue with the handling of strings in xarray.
It affects the coordinates with strings as labels such as :code:`radiating_dof` and :code:`influenced_dof`.
They can be stored in xarray either as NetCDF string objects, which can be written in a NetCDF file, or as Python strings stored as generic Python objects, which cannot be written in a NetCDF file.
The issue is that the xarray library sometimes changes from one to the other without warnings.
It leads to the error :code:`ValueError: unsupported dtype for netCDF4 variable: object` when trying to export a dataset.

This can be fixed by explicitely converting the strings to the right format when exporting the dataset::

    dataset.to_netcdf("dataset.nc",
                      encoding={'radiating_dof': {'dtype': 'U'},
                                'influenced_dof': {'dtype': 'U'}})

See also `this Github issue <https://github.com/mancellin/capytaine/issues/2>`_.

Saving the data as legacy Tecplot files
---------------------------------------

.. warning:: This feature is experimental.

The following code will write files named :code:`RadiationCoefficients.tec` and :code:`ExcitationForce.tec` in a format matching the one of Nemoh 2.0::

	from capytaine.io.legacy import write_dataset_as_tecplot_files
	write_dataset_as_tecplot_files("path/to/directory", dataset)

