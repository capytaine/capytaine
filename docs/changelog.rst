=========
Changelog
=========

------------------
New in version 0.5
------------------

Major changes
-------------

* Integration of a subset of Meshmagick into Capytaine as the :code:`mesh` submodule.
  Meshmagick is not a dependancy any more.
* Reorganization of some submodules::

  capytaine.mesh_collection -> capytaine.mesh.meshes_collection
  capytaine.symmetries -> capytaine.mesh.symmetries
  capytaine.cli -> capytaine.ui.cli
  capytaine.tools.vtk -> capytaine.ui.vtk
  capytaine.tools.mpl_free_surface_animation -> capytaine.ui.mpl_free_surface_animation
  capytaine.tools.import_export -> capytaine.io.legacy
  capytaine.tools.bemio -> capytaine.io.bemio
  meshmagick.mmio -> capytaine.io.mesh_loaders and capytaine.io.mesh_writers

* Reorganization of the internals of the solver :code:`Nemoh.py`. The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.
* :code:`CollectionOfMeshes` and its subclasses are now subclasses from :code:`tuple`.

Minor changes
-------------

* New function :code:`write_dataset_as_tecplot_files` in :code:`capytaine.tools`.
* New method :code:`tree_view()` for meshes to display the structure of hierarchical collections of meshes.
* Some body properties are stored in xarray dataset if they are available.
* Various improvements in :code:`geometric_bodies` submodule.

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

