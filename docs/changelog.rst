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
  The modules :code:`meshes_collection.py` and :code:`symmetries.py` have been moved to the :code:`mesh` submodule.

Minor changes
-------------

* Reorganization of the internals of the solver :code:`Nemoh.py`. The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.
* New function :code:`write_dataset_as_tecplot_files` in :code:`capytaine.tools`.

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

