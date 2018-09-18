=========
Changelog
=========

------------------
New in version 0.5
------------------

Major changes
-------------

* Experimental OpenMP parallelization of the computation of the influence matrices.
  The parallelization in :code:`solve_all` has been removed.

* Integration of a refactored subset of Meshmagick into Capytaine as the :code:`mesh` submodule.
  Meshmagick is not a dependancy any more.

* Reorganization of the submodules:

::

  capytaine.mesh_collection                  -> capytaine.mesh.meshes_collection
  capytaine.symmetries                       -> capytaine.mesh.symmetries
  capytaine.cli                              -> capytaine.ui.cli
  capytaine.tools.vtk                        -> capytaine.ui.vtk
  capytaine.tools.mpl_free_surface_animation -> capytaine.ui.mpl_free_surface_animation
  capytaine.tools.import_export              -> capytaine.io.legacy
  capytaine.tools.bemio                      -> capytaine.io.bemio
  meshmagick.mmio                            -> capytaine.io.mesh_loaders and capytaine.io.mesh_writers

Minor changes
-------------

Solver
~~~~~~

* Reorganization of the internals of the solver :code:`Nemoh.py` and :code:`NemohCore`.
  The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.
  The former has been replaced by a :code:`matrix_cache_size` optional argument (default value: 1).

* Support of :math:`\omega=0` and :math:`\omega=\infty` in the infinite depth case.

* The wavenumber is not computed in Fortran anymore.

Outputs
~~~~~~~

* Some body properties are stored in xarray dataset if they are available.
  New functions :code:`add_wavenumber_coords` and :code:`kochin_data_array` allow the storage of wavenumbers and Kochin function in the dataset.

* New functions :code:`separate_complex_values` and :code:`merge_complex_values`
  in :code:`capytaine.io.xarray` to better handle complex values when saving
  datasets.

* New function :code:`problems_from_dataset` to generate a list of problems from the coordinates of
  a xarray dataset.
  New method :code:`fill_dataset` in :code:`capytaine.Nemoh.Nemoh` using the above.

* New function :code:`write_dataset_as_tecplot_files()` in :code:`capytaine.tools` for legacy Tecplot output.

Meshes
~~~~~~

* Refactoring of the transformation methods (:code:`translate`, :code:`rotate`, :code:`mirror`, ...).

  * They are still in place by default, although they now return a reference to the modified object.
  * They can return a new object by passing the argument :code:`inplace=False` or by using the
    variants :code:`translated`, :code:`rotated`, :code:`mirrored`.
  * :code:`rotate` and :code:`rotated` requires an :code:`Axis` object as argument. Old behavior
    can be found in :code:`rotate_angles` and :code:`rotated_angles`.
  * :code:`get_immersed_part` is inplace by default. Use :code:`inplace=False` to return a new
    object.

* :code:`add_rotation_dof` now requires an Axis object.

* New method :code:`tree_view()` for meshes to display the structure of hierarchical collections of meshes.

* :code:`CollectionOfMeshes` and :code:`SymmetriBodies` are now subclasses from :code:`tuple`.
  New methods :code:`join_meshes` to merge several symmetric bodies with the same symmetries as a
  single symmetric body.

* Various improvements in :code:`geometric_bodies` submodule, especially for :code:`Rectangle` and :code:`RectangularParallelepiped`.
  They can now be generated with reflections symmetries instead of translation symmetries.
  New :code:`VerticalCylinder` class.

* Refactored mesh objects can be checked for equality and are hashable.
  The method is experimental and can be improved.

------------------
New in version 0.4
------------------

New features
------------

* Documentation and new usage examples.
* Computation of Kochin coefficients.
* Cleverer helper functions to define degrees of freedom.

Major changes
-------------

* Backward-incompatible change of the way the degrees of freedom are stored.

Minor changes
-------------

* Double precision computations.
* Improvement of :code:`assemble_dataset` for parametric studies.
* Support clipping of collections of meshes.
* Fixes in geometrical bodies generation.

