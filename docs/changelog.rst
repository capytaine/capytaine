=========
Changelog
=========

.. contents::
   :local:
   :depth: 1
   :backlinks: none

---------------------------------
New in version 1.4.2 (2022--)
---------------------------------

Bug fixes
~~~~~~~~~

* Raise error message when calling :meth:`~capytaine.bodies.bodies.FloatingBody.compute_hydrostatics()` without a center of mass defined (:pull:`207`).

* Fix bug when cropping body with a dof defined manually as a list of tuples (:issue:`204` and :pull:`206`).

Internals
~~~~~~~~~

* Replace the Fortran core by a git submodule pointing to `libDelhommeau <https://github.com/capytaine/capytaine/>`_ (:pull:`208`).
  Future developments of the Green function will take place there.

---------------------------------
New in version 1.4.1 (2022-09-05)
---------------------------------

Bug fixes
~~~~~~~~~

* Fix bug in hydrostatics of rigid bodies: the hydrostatic matrices were always assuming that the rotation dofs were defined around the :math:`(0, 0, 0)` point.
  The stiffness and inertia matrix are now invariant by horizontal translation of the body, as they should be. (:issue:`178` and :pull:`196`).

* Removed outdated volume/area methods from pre-defined bodies (:pull:`183`).

* Added symmetric realization and reflection to gdf mesh import (:issue:`186` and :pull:`187`).

* Fix some automatic mesh names (:pull:`195`)

* Fix ordering of the dofs when using :meth:`~capytaine.bodies.bodies.FloatingBody.assemble_regular_array()` (:issue:`198` and :pull:`199`)

* Return more explicit error message when the center of mass is missing for the computation of rigid-body hydrostatics (:pull:`201`).

* Return error message when trying to animate a body with a dof that has not been defined. Previously, undefined dofs were silently ignored. (:pull:`202`)


-------------------------------
New in version 1.4 (2022-07-07)
-------------------------------

Major changes
~~~~~~~~~~~~~

* The function that used to be called :code:`impedance` is now named :func:`~capytaine.post_pro.impedance.rao_transfer_function`.
  The new function :func:`~capytaine.post_pro.impedance.impedance` is the actual impedance matrix (:pull:`142`, :issue:`147`, :pull:`149`).

* The mass matrix of a floating body used to be denoted :code:`mass`. It is now denote :code:`inertia_matrix`.
  The attribute :code:`body.mass` is now used instead for the (scalar) mass of the body. (:pull:`165`)

* Implementation of :class:`~capytaine.bodies.predefined.spheres.Sphere` has changed.
  The use of symmetry is now controlled by the :code:`axial_symmetry` keyword argument.
  The :code:`clever` keyword argument is deprecated for :code:`Sphere` and should be replaced by the more explicit keyword arguments :code:`axial_symmetry`.
  Meanwhile, a bug has been fixed with its :code:`geometric_center` (:pull:`150`).

* The default linear solver is the direct solver and not the iterative solver GMRES, because it is more robust and more predictable.
  Nothing changes when users explicitely choose a linear solver. (:pull:`171`)

Bug fixes
~~~~~~~~~

* Fix major bug in impedance matrix and RAO computation: the sign of the dissipation matrix was wrong in previous versions (:issue:`102` and :pull:`140`).

* Fix major inaccuracy for deep panels or high frequencies, that is panels deeper than :math:`1.2\lambda` below the free surface where :math:`\lambda` is the wavelength (:issue:`38` and :pull:`156`)

* Wave directions in :code:`Nemoh.cal` are interpreted as degrees as they should be (and then converted to radians to be handled by the rest of the code). (:pull:`141`)

* Fix bug in rotations around axis that does not pass by (0, 0, 0) (:issue:`151` and :pull:`152`).

* New implementation of the mesh importer for :code:`hst` files. (:pull:`90`)
  It should be more robust and support more variants of the :code:`hst` mesh file format.

* Support for quadratures from `quadpy <https://pypi.org/project/quadpy/>`_ has been updated to support the version 0.16.16 of quadpy (:pull:`164`).

New features
~~~~~~~~~~~~

* Add method to compute some of the hydrostatic parameters such as volume, buoyancy center, wet surface area, hydrostatic stiffness, inertia matrix etc.
  :code:`compute_hydrostatics` method is created to return all hydrostatic parameters similar to :code:`meshmagick.hydrostatics.compute_hydrostatics` (:pull:`106`).
  By default, the hydrostatics are computed assuming a neutrally buoyant body (its mass is the displaced mass of water).
  Non-neutrally buoyant are partially supported, by setting the :code:`mass` attribute of the body (:pull:`166`)

* Add new parallelization using the `joblib <https://joblib.readthedocs.io>`_ library as a new optional dependency.
  The optional keyword-argument :code:`n_jobs` in the :meth:`~capytaine.bem.solver.BEMSolver.solve_all` and :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` controls the number of processes running in parallel (:pull:`136`). By default, this parallelisation is disabled (:pull:`172`).

* Refactor Delhommeau's method for the Green function evaluation. The size of the tabulation is not hard-coded anymore and can be changed by users. (:issue:`20` and :pull:`157`)

* Method :code:`show_matplotlib` can now colour mesh faces based on a specified scalar field (e.g. pressure) (:pull:`122`).

* The functions :func:`~capytaine.io.xarray.problems_from_dataset` and :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` accept a body alone as input.
  That is, one can use :code:`fill_dataset(test_matrix, body)` and not only :code:`fill_dataset(test_matrix, [body])` (:pull:`144`).

Documentation and error handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Improve feedback to users:
    * Print a warning if the user provides a :code:`wave_direction` that is not in the range [-2π, 2π]. (:pull:`141`)
    * Raise an error when the Green function evaluation returns a :code:`NaN` (:pull:`143`).
    * Improve message when the GMRES did not converge (:pull:`143`).
    * Raise an error when a body with an empty mesh is given to :code:`LinearPotentialFlowProblem` (:issue:`128` and :pull:`145`).
    * Print a warning when a key is unknown in the test matrix provided to :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` (:pull:`155`).
    * Raise an error if neither :code:`radiating_dof` (for radiation problems) nor :code:`wave_direction` (for diffraction problems) is provided in the test matrix in :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` (:pull:`155`).

* A new example using Haskind's relation has been added to the cookbook (:pull:`129`).

* Miscellaneous improvements of the documentation.


-------------------------------
New in version 1.3 (2021-10-07)
-------------------------------

Major changes
~~~~~~~~~~~~~

* The mesh are always "healed" when a new :code:`FloatingBody` is initialised
  (i.e. unused vertices are removed, degenerate triangles are removed, etc.).
  See for instance :issue:`46`.

* Implementation of symmetries for :code:`HorizontalCylinder` has changed.
  The cylinder is now a reflection of two halves containing translational
  repetition of half rings, instead of a translational repetition of symmetric
  ring. By default, only reflection symmetry is used. (:pull:`91`)
  The use of symmetries can be controlled with :code:`translation_symmetry` and
  :code:`reflection_symmetry` optional keyword arguments.
  The :code:`clever` keyword argument is deprecated for :code:`HorizontalCylinder`
  and should be replaced by the new more explicit keyword arguments above.


New features
~~~~~~~~~~~~

* Add method :code:`FloatingBody.from_meshio` to import `meshio <https://pypi.org/project/meshio/>`_ and `pygmsh <https://pypi.org/project/pygmsh/>`_ mesh objects (:pull:`62`)

* Add method :code:`FloatingBody.assemble_arbitrary_array` to make an array of bodies with arbitrary layout (:pull:`71`).

* Break out impedance from RAO to separate function (:issue:`61` and :pull:`63`).

* Method `problems_from_dataset` can now use a list of gravitational acceleration `g` values in the test matrix (:pull:`86`).

* Add example in cookbook for computing hydrostatics and mass properties with Meshmagick 2 (:pull:`70`).

Bug fixes
~~~~~~~~~

* Fix bug in free surface elevation computation when the number of faces in the free surface mesh is not a multiple of the chunk size, that is by default a multiple of 50 (:pull:`82`).

* The function :code:`assemble_dataset` did not support well the problems without a free surface. In the new version, such problems are explicitly ignored and a warning message is displayed. (:issue:`88` and :pull:`89`).

* Fix bug in some of the mesh readers/writers when using pathlib path objects (:pull:`87`).

* Function :code:`load_GDF` has been rewritten to accept any GDF file format (:pull:`97`).

Internal and development
~~~~~~~~~~~~~~~~~~~~~~~~

* Easier installation of optional dependencies via :code:`pip install -e .[extra]` and :code:`pip install -e .[develop]` (:pull:`96`).

* Use pytest skipif to skip tests if optional dependencies are not installed (:pull:`68`).

---------------------------------
New in version 1.2.1 (2021-04-14)
---------------------------------

* Minor bug fixes,
  including :issue:`37`
  and :issue:`56` (thanks to Ryan Coe).

* Add a warning when a panel is on the free surface
  (see :issue:`29` and :issue:`50`)

-------------------------------
New in version 1.2 (2020-04-24)
-------------------------------

* Experimental implementation of higher order quadratures for the integration of
  the Green function on the mesh. Default behavior is still the first order
  integration as in Nemoh.

* Add method :code:`FloatingBody.animate` to quickly visualize the motion of a body
  and add method :code:`Animation.embed_in_notebook` to embed animations in Jupyter
  notebooks.

* Keep the order of the dofs in `xarray`'s Datasets.
  This patch uses the CategoricalIndex feature of `xarray` which was buggy
  before version 0.15.1 of `xarray`. Thus this minimal version is now required.

* Add missing Kochin function for the diffraction.
  (See :issue:`22`.)
  In previous version the variable named :code:`kochin` in the dataset was only the
  Kochin function for the radiated waves. A new variable names
  :code:`kochin_diffraction` has been added. The existing variable :code:`kochin` has not
  been renamed, for backward compatibility, but might be in a future release of
  Capytaine.

* Improvement of caching to limit RAM usage for large problems.

* Make optional the dependency to graphical packages (`matplotlib` and `vtk`).
  They were causing issues to some users.

* :code:`problems_and_results.py` has been rewritten to be slightly more readable and
  remove the dependency to `attrs`.

-------------------------------
New in version 1.1 (2019-09-24)
-------------------------------

Major changes
~~~~~~~~~~~~~

* Refactoring of the implementation of the solver.
  The new implementation separates the solver itself from the evaluation of the
  Green function and the matrix building engine.
  This more modular structure allows user to choose separately the Green
  function and the matrix engine that they want to use.

  The former API (:code:`Nemoh()` object) has been kept for backward compatibility.
  In most cases, replacing :code:`Nemoh()` by :code:`BEMSolver()` is sufficient
  to migrate to the new structure.

  See :doc:`user_manual/resolution` for the full documentation of the new structure
  and :doc:`user_manual/cookbook` for examples.


* Add Xie's variant of Delhommeau's Green function
  :class:`~capytaine.green_functions.delhommeau.XieDelhommeau` [X18]_.


* The option `cache_rankine_matrices` has been removed. It was impeding the
  performance and modularity of the code for a very low gain. It might be
  reimplemented in a future version if there is really a need for it.

Minor changes
~~~~~~~~~~~~~

* Minor performance improvements.

* Fix Github issue #18.

* Improve test suite.

---------------------------------
New in version 1.0.1 (2019-03-28)
---------------------------------

Minor changes
~~~~~~~~~~~~~

* Fix compilation flags for OpenMP

* Minor corrections in the documentation.

-------------------------------
New in version 1.0 (2019-03-14)
-------------------------------

Major changes
~~~~~~~~~~~~~

* The :code:`angle` parameter has been renamed to the more accurate name
  :code:`wave_direction`.

* Most of the modules have been reorganized in several packages. See the
  :doc:`developer_manual/overview` for some details.

* Test compatibility of the code with Python 3.7 and numpy 1.16.

* Remove a couple of unmaintained or unfinished submodules.

Minor changes
-------------

General
~~~~~~~

* Many improvements of the documentation.

* Reorganization of some of the tests.

* Various small performance improvement.

Mesh and bodies
~~~~~~~~~~~~~~~

* Rename :code:`center` into either :code:`geometric_center` or
  :code:`center_of_mass` depending of the case.

* New method for geometric bodies :code:`rotate_around_center_to_align_vectors`
  replacing :code:`rotate_to_align_axes`.

* Add methods :code:`sliced_by_plane` and :code:`minced` for hierarchical
  decomposition of the mesh.

* Symmetric meshes classes have been renamed::

    ReflectionSymmetry -> ReflectionSymmetricMesh
    etc.

* Plane are now oriented: they are equal only if their normal point in the same
  direction.

Solver
~~~~~~

* Store solver settings in output dataset.

* Rename setting :code:`use_symmetries` into :code:`hierarchical_toeplitz_matrices`.

* Fix bugs and improve implementation of the Adaptive Cross Approximation.

-------------------------------
New in version 0.6 (2019-02-11)
-------------------------------

Major changes
~~~~~~~~~~~~~

* Full rewrite of the matrices and linear solvers implementation.
  All the relevant code is now in the submodule :code:`capytaine.matrices`.

* Refactored implementation of block Toeplitz matrices, block symmetric Toeplitz
  matrices and block circulant matrices.
  Refactoring of the block diagonalization of block circulant matrices through
  FFT.

* Low rank approximation of the matrices with Adaptive Cross Approximation for
  the use of hierarchical matrices.

* Option to solve the linear system with GMRES instead of a direct solver.

* Refactoring of the 3D animation module for animation of the body motions,
  animated colormap of the pressure, free-surface elevation and export as a
  video. See cookbook for an example of the new API.

Minor changes
~~~~~~~~~~~~~

General
-------

* Reorganization of the :code:`pytest` directory.

* Add an experimental :code:`capytaine.tools.rao` module to compute Response Amplitude
  Operators.

* Various bug fixes and improvements of the documentation.

Solver
------

* More options to set the behavior of the solver at run time :code:`Nemoh` (use
  of symmetries, use of caching, choice of the linear solver, ...).
  See its docstring for details.

* Change of default behavior: the solver stores the details in the :code:`Result`
  container when using :code:`solve`, not when using :code:`solve_all` or
  :code:`fill_dataset`.

* The water density can be specified in the test matrix when using
  :code:`fill_dataset`.

* Function :code:`kochin_dataset` to build a xarray of Kochin function.

* Add the option :code:`chunk_size` to the computation of the free surface
  elevation in order to limit the RAM consumption.

* Minor refactoring of the solver and the computation of the Green function.

Meshes and bodies
-----------------

* CollectionOfMeshes is not a subclass of Tuple anymore.

* New method :code:`assemble_regular_array` to build an array of identical bodies.

* Harmonize naming of functions that are not in-place: :code:`symmetrize -> symmetrized`, :code:`merge -> merged`

* Refactoring of the internals of the mesh clipper. New :code:`clip` and :code:`clipped` methods for meshes and bodies.
  When a body is clipped with :code:`clip` or :code:`keep_immersed_part`, the dofs are updated.

* Change naming of dof when bodies are joined: :code:`body_name__dof_name` instead of :code:`body_name_dof_name`.

* The combination of bodies with :code:`+` is associative with respect to the
  names of the dofs.

* Minor improvements of meshes and bodies :code:`repr`.

---------------------------------
New in version 0.5.1 (2018-10-17)
---------------------------------

* Minor bugs fixes.

* Small performance improvements.

* Update documentation.

-------------------------------
New in version 0.5 (2018-09-18)
-------------------------------

Major changes
~~~~~~~~~~~~~

* Experimental OpenMP parallelization of the computation of the influence matrices.
  The parallelization in :code:`solve_all` has been removed.

* Integration of a refactored subset of Meshmagick into Capytaine as the :code:`mesh` submodule.
  Meshmagick is not a dependency any more.

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
~~~~~~~~~~~~~

Solver
------

* Reorganization of the internals of the solver :code:`Nemoh.py` and :code:`NemohCore`.
  The initialization options :code:`keep_matrices` and :code:`max_stored_exponential_decompositions` have been removed.
  The former has been replaced by a :code:`matrix_cache_size` optional argument (default value: 1).

* Support of :math:`\omega=0` and :math:`\omega=\infty` in the infinite depth case.

* The wavenumber is not computed in Fortran anymore.

Outputs
-------

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
------

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

* :code:`CollectionOfMeshes` and :code:`SymmetricBodies` are now subclasses from :code:`tuple`.
  New methods :code:`join_meshes` to merge several symmetric bodies with the same symmetries as a
  single symmetric body.

* Various improvements in :code:`geometric_bodies` submodule, especially for :code:`Rectangle` and :code:`RectangularParallelepiped`.
  They can now be generated with reflections symmetries instead of translation symmetries.
  New :code:`VerticalCylinder` class.

* Refactored mesh objects can be checked for equality and are hashable.
  The method is experimental and can be improved.

-------------------------------
New in version 0.4 (2018-08-04)
-------------------------------

New features
~~~~~~~~~~~~

* Documentation and new usage examples.
* Computation of Kochin coefficients.
* Cleverer helper functions to define degrees of freedom.

Major changes
~~~~~~~~~~~~~

* Backward-incompatible change of the way the degrees of freedom are stored.

Minor changes
~~~~~~~~~~~~~

* Double precision computations.
* Improvement of :code:`assemble_dataset` for parametric studies.
* Support clipping of collections of meshes.
* Fixes in geometrical bodies generation.
