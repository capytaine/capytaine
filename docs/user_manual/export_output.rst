==============
Export outputs
==============


The ``Dataset`` output format is a standard object from Xarray and all methods from `Xarray's manual <https://docs.xarray.dev/en/stable/user-guide/io.html>`_ can be used to manipulate and export it.
On top of that, Capytaine provides some wrappers functions to simplify the export into a NetCDF file, as well as into other formats more specific to hydrodynamics.

.. contents:: Contents

Exporting hydrodynamic dataset
------------------------------

NetCDF format
~~~~~~~~~~~~~

The xarray dataset produced by :func:`assemble_dataset <capytaine.io.xarray.assemble_dataset>` (or :meth:`fill_dataset <capytaine.bem.solver.BEMSolver.fill_dataset>`) has a structure close to the NetCDF file format and can easily be saved to this format by :func:`~capytaine.io.xarray.export_dataset`::

    cpt.export_dataset("path/to/dataset.nc", dataset, format="netcdf")


.. note::
    The netCDF standard does not handle **complex numbers**.

    Capytaine is using a non-standard representation of complex numbers by
    transforming all complex-valued arrays of shape ``(...)`` into real-valued
    arrays of shape ``(..., 2)``` using the functions
    :func:`~capytaine.io.xarray.separate_complex_values` and
    :func:`~capytaine.io.xarray.merge_complex_values` (done automatically by
    :func:`~capytaine.io.xarray.export_dataset`).


    See also https://github.com/PlasmaFAIR/nc-complex for more context and alternatives.

.. note::
    Exporting more outputs such as pressure field on the hull in a NetCDF
    file is considered in the future.
    See https://github.com/capytaine/capytaine/issues/520 for examples of such outputs.

Wamit format
~~~~~~~~~~~~

The hydrodynamic results from a Capytaine ``xarray.Dataset`` can be exported into WAMIT-compatible text files (``.1``, ``.3``, ``.3fk``, ``.3sc``, ``.hst``) using::

    cpt.export_dataset("problem_name", dataset, format="wamit", exports=("1", "3", "3fk", "3sc", "hst"))

This will produce the following files (depending on the fields present in the dataset and the flags passed to the optional ``exports`` argument):

* ``problem_name.1`` for added mass and radiation damping coefficients,

* ``problem_name.3`` for total excitation forces (Froude-Krylov + diffraction),

* ``problem_name.3fk`` for Froude-Krylov forces only,

* ``problem_name.3sc`` for diffraction forces only.

* ``problem_name.hst`` for hydrostatics results (if supported)

Invalid or unavailable exports are skipped with a warning.

The length scale used for normalization in WAMIT data is taken by default as :math:`1` meter.

.. note::
    These exports require that the ``forward_speed`` in the dataset is zero.
    If not, a ``ValueError`` is raised to avoid exporting inconsistent results.


Nemoh format
~~~~~~~~~~~~

The following code will write files named :code:`RadiationCoefficients.tec` and :code:`ExcitationForce.tec` in a format roughly matching the one of Nemoh 2::

    cpt.export_dataset("path/to/result_dir/", dataset, format="nemoh")

This feature is still experimental. Please report issues encountered with this.


Excel format
~~~~~~~~~~~~

Export to Excel format is not currently built in :func:`~capytaine.io.xarray.export_dataset`.
This section is meant to show an example of exporting to a format that is not explicitly implemented in Capytaine.
We use here the ``openpyxl`` library (that can be installed with ``pip install openpyxl``) to export a dataset to Excel format::

    dataset[["added_mass", "radiation_damping"]].to_dataframe().to_excel("radiation_data.xlsx")

    from capytaine.io.xarray import separate_complex_values
    separate_complex_values(dataset[["Froude_Krylov_force", "diffraction_force"]]).to_dataframe().to_excel("diffraction_data.xlsx")

For convenience, the radiation and diffraction data have been stored in separate files.
Since this export method poorly supports complex number, the :func:`separate_complex_values <capytaine.io.xarray.separate_complex_values>` has been used to transform them to a pair of real numbers, as discussed for NetCDF export above.


Saving the rotation center of rigid bodies in NetCDF files
----------------------------------------------------------

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
  cpt.export_dataset("dataset.nc", dataset, format="netcdf")

The support for this in Capytaine should be improved in the future.


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
