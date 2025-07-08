==============
Export outputs
==============

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

    import xarray as xr
    from capytaine.io.xarray import merge_complex_values
    dataset = merge_complex_values(xr.open_dataset("path/to/dataset.nc"))


String format
~~~~~~~~~~~~~

There is an issue with the handling of strings in xarray.
It affects the coordinates with strings as labels such as :code:`radiating_dof` and :code:`influenced_dof`.
They can be stored in xarray either as NetCDF string objects, which can be written in a NetCDF file, or as Python strings stored as generic Python objects, which cannot be written in a NetCDF file.
The issue is that the xarray library sometimes changes from one to the other without warnings.
It leads to the error :code:`ValueError: unsupported dtype for netCDF4 variable: object` when trying to export a dataset.

This can be fixed by explicitly converting the strings to the right format when exporting the dataset::

    separate_complex_values(dataset).to_netcdf(
      "dataset.nc",
      encoding={'radiating_dof': {'dtype': 'U'},
                'influenced_dof': {'dtype': 'U'}}
    )

See also `this Github issue <https://github.com/capytaine/capytaine/issues/2>`_.


Saving the rotation center of rigid bodies
------------------------------------------

Saving rotation hydrodynamic coefficients without explicitly defining the rotation axes can be ambiguous and can lead to confusion downstream.
While this is not done automatically by Capytaine at the moment, it can be added to the dataset manually.
The example below, which is an extension of the :doc:`quickstart` example, saves the rotation centers of a multibody problem in a way that is understood notably by `BEMRosetta <https://github.com/BEMRosetta/BEMRosetta>`_::

  import numpy as np
  import xarray as xr
  import capytaine as cpt

  body_1 = cpt.FloatingBody(
              mesh=cpt.mesh_sphere(center=(0, 0, 0)),
              dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
              center_of_mass=(0, 0, 0),
              name="my_sphere",
          )
  body_1.inertia_matrix = body_1.compute_rigid_body_inertia()
  body_1.hydrostatic_stiffness = body_1.immersed_part().compute_hydrostatic_stiffness()
  # If you have several rigid bodies, copy the code above to define "body_2", "body_3", etc.

  list_of_bodies = [body_1]  # Replace "[body_1]" by "[body_1, body_2, body_3]" for multibody problem.

  all_bodies = cpt.FloatingBody.join_bodies(*list_of_bodies).immersed_part()

  # Set up parameters
  test_matrix = xr.Dataset({
          "omega": np.linspace(0.1, 2.0, 20),  # Can also specify "period", "wavelength" or "wavenumber"
          "wave_direction": np.linspace(0, np.pi, 3),
          "radiating_dof": list(all_bodies.dofs),
          })

  # Do the resolution
  solver = cpt.BEMSolver()
  dataset = solver.fill_dataset(test_matrix, all_bodies)

  dataset.coords["rigid_body_component"] = [body.name for body in list_of_bodies]
  dataset["rotation_center"] = (["rigid_body_component", "point_coordinates"], [body.rotation_center for body in list_of_bodies])
  dataset["center_of_mass"] = (["rigid_body_component", "point_coordinates"], [body.center_of_mass for body in list_of_bodies])

  # Export to NetCDF file
  from capytaine.io.xarray import separate_complex_values
  separate_complex_values(dataset).to_netcdf("dataset.nc",
                                             encoding={'radiating_dof': {'dtype': 'U'},
                                                       'influenced_dof': {'dtype': 'U'}})

The support for this in Capytaine should be improved in the future.

Exporting to Excel
------------------

The example below uses the ``openpyxl`` library (that can be installed with ``pip install openpyxl``) to export a dataset to Excel format::

    dataset[["added_mass", "radiation_damping"]].to_dataframe().to_excel("radiation_data.xlsx")

    from capytaine.io.xarray import separate_complex_values
    separate_complex_values(dataset[["Froude_Krylov_force", "diffraction_force"]]).to_dataframe().to_excel("diffraction_data.xlsx")

For convenience, the radiation and diffraction data have been stored in separate files.
Since this export method poorly supports complex number, the :func:`separate_complex_values <capytaine.io.xarray.separate_complex_values>` has been used to transform them to a pair of real numbers, as discussed for NetCDF export above.


Saving the hydrostatics data of rigid body(ies) in Nemoh's format
-----------------------------------------------------------------

For a rigid body, or a set of several rigid bodies, the following information can be saved as written by Nemoh's and read by BEMIO to produce :code:`.h5` files for WEC-Sim:

- Hydrostatic stiffness matrix,
- Centre of gravity,
- Centre of buoyancy,
- Displacement volume

They are stored in two files (:code:`Hydrostatics.dat` and :code:`KH.dat`) for each body, using the following syntax::

    from capytaine.io.legacy import export_hydrostatics
    export_hydrostatics("directory_to_save_hydrostatics_data", body)

for a single rigid body or, e.g.,::

    from capytaine.io.legacy import export_hydrostatics
    export_hydrostatics("directory_to_save_hydrostatics_data", [body_1, body_2, body_3])

for several rigid bodies.

In order to use this function, please ensure that the body's centre of gravity has been defined correctly and the following methods have been called on the :code:`FloatingBody` object before passing it to :func:`export_hydrostatics <capytaine.io.legacy.export_hydrostatics>`::

  body.add_all_rigid_body_dofs()
  body.inertia_matrix = body.compute_rigid_body_inertia()
  body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()


Saving the data as legacy Tecplot files
---------------------------------------

.. warning:: This feature is experimental.

The following code will write files named :code:`RadiationCoefficients.tec` and :code:`ExcitationForce.tec` in a format matching the one of Nemoh 2.0::

	from capytaine.io.legacy import write_dataset_as_tecplot_files
	write_dataset_as_tecplot_files("path/to/directory", dataset)


Exporting to WAMIT format
-------------------------

The hydrodynamic results from a Capytaine ``xarray.Dataset`` can be exported into WAMIT-compatible text files (``.1``, ``.3``, ``.3fk``, ``.3sc``, ``.hst``) using::

    from capytaine.io.wamit import export_to_wamit
    export_to_wamit(dataset, "problem_name", exports=("1", "3", "3fk", "3sc", "hst"))

This will produce the following files (depending on the fields present in the dataset):

* ``problem_name.1`` for added mass and radiation damping coefficients,

* ``problem_name.3`` for total excitation forces (Froude-Krylov + diffraction),

* ``problem_name.3fk`` for Froude-Krylov forces only,

* ``problem_name.3sc`` for diffraction forces only.

* ``problem_name.hst`` for hydrostatics results (if supported)

Invalid or unavailable exports are skipped with a warning.

The length scale used for normalization in WAMIT data is taken by default as
:math:`1` meter.

.. note::
    These exports require that the ``forward_speed`` in the dataset is zero.
    If not, a ``ValueError`` is raised to avoid exporting inconsistent results.


Building a dataset from Bemio
-----------------------------

A DataFrame or a Dataset can also be created from data structures generated
using the `Bemio <https://wec-sim.github.io/bemio/>`_ package, which reads
hydrodynamic output data from NEMOH, WAMIT, and AQWA. This allows for Capytaine
post-processing of hydrodynamic data generated from other BEM codes.

Bemio does not come packaged with Capytaine and needs to to be installed independently.
Note that `the base repository of Bemio <https://github.com/WEC-Sim/bemio/>`_ has been
archived and is only compatible with Python 2.7.x, so using a Python 3 compatible fork is
recommended, available `here <https://github.com/mancellin/bemio>`_ or installed with::

  pip install git+https://github.com/mancellin/bemio.git

To build the xarray dataset using Capytaine, the output files from the BEM program in
question must be read into a Bemio :code:`data_structures.ben.HydrodynamicData` class, which is
then called by `assemble_dataframe` or `assemble_dataset`. For example, to
create an xarray dataset from a WAMIT :code:`.out` file::

  from bemio.io.wamit import read as read_wamit
  import capytaine as cpt
  bemio_data = read_wamit("myfile.out")
  my_dataset = cpt.assemble_dataset(bemio_data, hydrostatics=False)

.. warning:: The created dataset will only contain quantities that can be directly calculated
             from the values given in the original dataset. Because of this, there may be minor
             differences between the variable names in an xarray dataset build with Bemio and one created
             using :code:`LinearPotentialFlowResult`, even though the format will be identical. For
             example, WAMIT :code:`.out` files do not contain the radii of gyration needed to calculate
             the moments of inertia, so the `my_dataset['inertia_matrix']` variable would not be included
             in the above example since the rigid body mass matrix cannot be calculated.
