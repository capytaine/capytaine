=========
Changelog
=========

------------------
New in version 0.5
------------------

Minor changes
-------------

* Reorganization of the internals of the solver. The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.

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

