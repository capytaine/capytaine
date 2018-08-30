=========
Changelog
=========

------------------
New in version 0.5
------------------

New features
------------

* New function :code:`write_dataset_as_tecplot_files()` in :code:`capytaine.tools`.

* New method :code:`tree_view()` for meshes to display the structure of hierarchical collections of meshes.

* Some body properties are stored in xarray dataset if they are available. New functions :code:`add_wavenumber_coords` and :code:`kochin_data_array` allow the storage of wavenumbers and Kochin function in the dataset.

* New function :code:`problems_from_dataset` to generate a list of problems from the coordinates of
  a xarray dataset.
  New method :code:`fill_dataset` in :code:`capytaine.Nemoh.Nemoh` using the above.

Major changes
-------------

* Integration of a refactored subset of Meshmagick into Capytaine as the :code:`mesh` submodule.
  Meshmagick is not a dependancy any more.
* Reorganization of some submodules:

::

  capytaine.mesh_collection                  -> capytaine.mesh.meshes_collection
  capytaine.symmetries                       -> capytaine.mesh.symmetries
  capytaine.cli                              -> capytaine.ui.cli
  capytaine.tools.vtk                        -> capytaine.ui.vtk
  capytaine.tools.mpl_free_surface_animation -> capytaine.ui.mpl_free_surface_animation
  capytaine.tools.import_export              -> capytaine.io.legacy
  capytaine.tools.bemio                      -> capytaine.io.bemio
  meshmagick.mmio                            -> capytaine.io.mesh_loaders and capytaine.io.mesh_writers

* Refactoring of the transformation methods (:code:`translate`, :code:`rotate`, :code:`mirror`, ...).

  * They are still in place by default, although they now return a reference to the modified object.
  * They can return a new object by passing the argument :code:`inplace=False` or by using the
    variants :code:`translated`, :code:`rotated`, :code:`mirrored`.
  * :code:`rotate` and :code:`rotated` requires an :code:`Axis` object as argument. Old behavior
    can be found in :code:`rotate_angles` and :code:`rotated_angles`.
  * :code:`get_immersed_part` is inplace by default. Use :code:`inplace=False` to return a new
    object.

* :code:`add_rotation_dof` now requires an Axis object.

* Reorganization of the internals of the solver :code:`Nemoh.py` and :code:`NemohCore`.
  The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.

* :code:`CollectionOfMeshes` and :code:`SymmetriBodies` are now subclasses from :code:`tuple`.

* Improvements in :code:`geometric_bodies` submodule, especially for :code:`Rectangle` and
  :code:`RectangularParallelepiped`.
  They can now be generated with reflections symmetries instead of translation symmetries.
  The argument :code:`clever` has been replaced by more relevantly named arguments.

* Refactored mesh objects can be checked for equality and are hashable. The method is experimental
  and can be improved.

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

