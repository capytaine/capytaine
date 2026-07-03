==============
Export outputs
==============


The ``Dataset`` output format is a standard object from Xarray and all methods from `Xarray's manual <https://docs.xarray.dev/en/stable/user-guide/io.html>`_ can be used to manipulate and export it.
On top of that, Capytaine provides some wrappers functions to simplify the export into a NetCDF file, as well as into other formats more specific to hydrodynamics.

.. contents:: Contents

NetCDF format
-------------

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
------------

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
------------

The results from a Capytaine ``xarray.Dataset`` can be exported into Nemoh-computible text files (:code:`Hydrostatics.dat`, :code:`KH.dat`, :code:`RadiationCoefficients.tec` and :code:`ExcitationForce.tec`) using::

    cpt.export_dataset("path/to/result_dir/", dataset, format="nemoh")

This feature is still experimental. Please report issues encountered with this.


Alternatively, :code:`Hydrostatics.dat` and :code:`KH.dat` can be exported with the following syntax::

    from capytaine.io.legacy import export_hydrostatics
    body.hydrostatic_stiffness = body.immersed_part().compute_hydrostatic_stiffness()
    export_hydrostatics("directory_to_save_hydrostatics_data", body)

for a single rigid body or, e.g.,::

    from capytaine.io.legacy import export_hydrostatics
    body_1.hydrostatic_stiffness = body_1.immersed_part().compute_hydrostatic_stiffness()
    body_2.hydrostatic_stiffness = body_2.immersed_part().compute_hydrostatic_stiffness()
    body_3.hydrostatic_stiffness = body_3.immersed_part().compute_hydrostatic_stiffness()
    export_hydrostatics("directory_to_save_hydrostatics_data", body_1 + body_2 + body_3)

for several rigid bodies.

Hydrostatics
~~~~~~~~~~~~

Unlike other formats, exporting the hydrostatics in legacy Nemoh format currently still requires a few extra steps.

For a rigid body, or a set of several rigid bodies, the following information can be saved as written by Nemoh's and read by BEMIO to produce :code:`.h5` files for WEC-Sim:

- Hydrostatic stiffness matrix,
- Centre of gravity,
- Centre of buoyancy,
- Displacement volume


In order to use this function, please ensure that the body's centre of gravity has been defined as well as the six rigid body degrees of freedom::

    body = cpt.FloatingBody(
            ...,
            dofs=cpt.rigid_body_dofs(rotation_center=...),
            center_of_mass=...,
            )


Excel format
------------

Export to Excel format is not currently built in :func:`~capytaine.io.xarray.export_dataset`.
This section is meant to show an example of exporting to a format that is not explicitly implemented in Capytaine.
We use here the ``openpyxl`` library (that can be installed with ``pip install openpyxl``) to export a dataset to Excel format::

    dataset[["added_mass", "radiation_damping"]].to_dataframe().to_excel("radiation_data.xlsx")

    from capytaine.io.xarray import separate_complex_values
    separate_complex_values(dataset[["Froude_Krylov_force", "diffraction_force"]]).to_dataframe().to_excel("diffraction_data.xlsx")

For convenience, the radiation and diffraction data have been stored in separate files.
Since this export method poorly supports complex number, the :func:`separate_complex_values <capytaine.io.xarray.separate_complex_values>` has been used to transform them to a pair of real numbers, as discussed for NetCDF export above.
