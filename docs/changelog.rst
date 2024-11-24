=========
Changelog
=========

.. contents::
   :local:
   :depth: 1
   :backlinks: none

---------------------------------
New in version 2.2.1 (2024-11-18)
---------------------------------

Minor change
~~~~~~~~~~~~

* More efficient implementation of the mesh connected-component clustering algorithm (:pull:`603`).

Bug fixes
~~~~~~~~~

* Lid meshes on the free surface do not cause errors when clipped.
  Also empty lid meshes are properly handled when initializing or clipping a mesh
  (:issue:`573` and :pull:`575`).

* GDF meshes are accepted in the alternate format now.
  Meshes files can list points in either 3 x 4*nPanels or a 12 x nPanels format.
  (:issue:`540` and :pull:`585`).

* When filling a test matrix with both diffraction problems and radiation
  problems, zero and infinite frequencies can now be provided. (Previously, the
  computation was failing because these frequencies are not defined for
  diffraction problems.) (:pull:`587`)

* Radiation damping at infinite frequency is now zero instead of infinity.
  When forward speed is non-zero, added mass and radiation dampings at zero encounter frequency are NaN.
  (:pull:`588`)

* User does not need to import ``pyplot`` themself before running `show_matplotlib()` (:pull:`592`)

* Fixes usage of ``ReflectionSymmetricMesh`` with direct solver (:issue:`593` and :pull:`594`).

* Do not recompute the same
  :meth:`~capytaine.bodies.bodies.FloatingBody.first_irregular_frequency_estimate``
  for the same body several times.
  Also better expose the ``_check_wavelength`` option to skip wavelength check,
  including irregular frequency estimation. (:issue:`601` and :pull:`602`).

* Fix bug in the :math:`r`-range of the tabulation of the Green function
  (:issue:`538` and :pull:`611`).

-------------------------------
New in version 2.2 (2024-07-08)
-------------------------------

Major changes
~~~~~~~~~~~~~

* **New feature: lid-based irregular frequencies removal**.
  Add ``lid_mesh`` argument to :class:`~capytaine.bodies.bodies.FloatingBody`
  for irregular frequencies removal (:pull:`521`).
  Add :meth:`~capytaine.meshes.meshes.Mesh.generate_lid` method to generate
  lids (:pull:`477`) and :meth:`~capytaine.meshes.meshes.Mesh.extract_lid`
  method to extract lids from exernally defined meshes (:pull:`559`).
  Add a warning to the user if irregular frequencies can be expected (:pull:`564`).

* The compiled Fortran extension is not split into a ``Delhommeau`` and a ``XieDelhommeau`` version anymore.
  The same effect is now achieved by the run-time parameter ``gf_singularities`` of the class :class:`~capytaine.green_functions.delhommeau.Delhommeau` (:pull:`475`).
  (The class :class:`~capytaine.green_functions.delhommeau.XieDelhommeau` is kept for backward compatibility.).
  The new default method in infinite depth is ``gf_singularities="low_freq"`` (formerly ``XieDelhommeau``) instead of ``gf_singularities="high_freq"``.
  The new one is expected to be more accurate near the surface and at low frequency (:pull:`566`)
  The finite depth Green function is always computed using the ``low_freq`` variant, so the ``gf_singularities`` parameter has no effect in finite depth. (:pull:`507`).
  The tabulation stores the data of both variants and is thus slightly longer to initialize and slightly larger to store in memory (:pull:`543`).

* Experimental support for panels on the free surface, when using ``gf_singularities="low_freq"``.  (:pull:`419`)

Minor changes
~~~~~~~~~~~~~

* Remove mesh resolution warning when the frequency is infinite (or the wavelength is zero) (:pull:`511`).

* When computing without a tabulation (``tabulation_nr=0`` or ``tabulation_nz=0``), the value of ``tabulation_nb_integration_points`` is actually used to compute Guével-Delhommeau exact formulation of the Green function. Previously, it was only used when precomputing a tabulation (:pull:`514`).

* Add a new variant of the Green function integration ``gf_singularities="low_freq_with_rankine_part"`` as an experimental more accurate version of the ``low_freq`` variant (:pull:`510`).

* Add a ``tabulation_cache_dir`` parameter to :class:`~capytaine.green_functions.delhommeau.Delhommeau` to choose the directory in which the tabulation is saved on disk. If ``None`` is provided instead, the tabulation is not saved on disk and is recomputed at each initialization of the class. Also, if this parameter is not set, look for the ``CAPYTAINE_CACHE_DIR`` environment variable and use it to save the tabulation if it exists. (:pull:`516`).

* Meshio objects can be directly passed to :func:`~capytaine.io.meshes_loaders.load_mesh` to get a Capytaine mesh (:pull:`555`).

* Load gmsh v4 format .msh file using :code:`cpt.load_mesh()` (when meshio is installed) (:pull:`556`)


Bug fixes
~~~~~~~~~

* Always use an odd number of points for integration with Simpson rule (:pull:`515`). This bug was partly responsible for some high-frequency inaccuracy (:issue:`298`).

* :func:`~capytaine.meshes.predefined.cylinders.mesh_vertical_cylinder` used to return only half of the mesh when called with ``reflection_symmetry=True`` (:issue:`529` and :pull:`530`).

* Providing the frequency as a scalar coordinate in the test matrix does not result in the value being ignored anymore (:issue:`547` and :pull:`548`).

* Improve exception message when giving an unknown ``radiating_dof`` to a :class:`~capytaine.bem.problems_and_results.RadiationProblem` (:pull:`549`).

* Fix issue due to breaking change in linear solver broadcasting in Numpy 2.0 (:issue:`550`).

* Remove warning mentioning missing divergence for rigid body dofs when computing hydrostatics (:pull:`487` and :pull:`570`)

Internals
~~~~~~~~~

* Update test environments used in noxfile and add ``editable_install_requirements.txt``. (:pull:`498`)

* Rename ``tabulation_method`` parameter of :class:`~capytaine.green_functions.delhommeau.Delhommeau` as the more descriptive ``tabulation_grid_shape``, and similarly for internal variables. (:pull:`503`)

* Add :func:`~capytaine.meshes.properties.connected_components` and :func:`~capytaine.meshes.properties.connected_components_of_waterline` to split a mesh into connected components. (:pull:`554`)

-------------------------------
New in version 2.1 (2024-04-08)
-------------------------------

Major changes
~~~~~~~~~~~~~

* **New feature: Approximate forward speed for single rigid body**.
  A ``forward_speed`` parameter can now be provided to :class:`~capytaine.bem.problems_and_results.LinearPotentialFlowProblem` (or to the test matrix when using :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset`) to compute the excitation force, added mass and radiation damping with forward speed of the body in the :math:`x` direction.
  Note that the :class:`~capytaine.bem.problems_and_results.RadiationProblem` now accept a ``wave_direction`` parameter, which is only used when ``forward_speed`` is non zero to compute the encounter frequency.
  See the theory manual for references. (:pull:`376`)

* Add `rich <https://rich.readthedocs.io>`_ as a dependency and improve formatting of the console output.
  Add :func:`~capytaine.ui.rich.set_logging` function to quickly set up logging with `rich`.
  :meth:`~capytaine.bem.solver.BEMSolver.solve_all` and :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` now display a progress bar (unless turn off by the ``progress_bar`` argument). (:pull:`382`)

* Reimplement computation of added mass and radiation damping in infinite depth with zero or infinite frequency. (:pull:`385` and :pull:`485`)
  When using forward speed, the added mass and radiation damping are undefined, but the forces can still be computed. (:pull:`483`)

* Implement direct method (source-and-dipole formulation) in obtaining velocity potential solutions.
  The direct method can be used instead of the default indirect method by setting the ``method`` argument of :meth:`~capytaine.bem.solver.BEMSolver.solve`, :meth:`~capytaine.bem.solver.BEMSolver.solve_all` or :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` (:pull:`420`)

* Add new shape for the grid used for the tabulation, based on the one used in Nemoh version 3.
  User can choose to use the Nemoh 3 grid shape (by default) or the former one by setting the ``tabulation_method`` parameter of :class:`~capytaine.green_functions.delhommeau.Delhommeau`.
  The new grid shape allows to set both the number of points (with ``tabulation_nr`` and ``tabulation_nz``) and the extent of the tabulation (with ``tabulation_rmax`` and ``tabulation_zmin``).
  The new default tabulation might lead to slightly different results, which are likely more accurate in the new version.
  (:pull:`439`)

Minor changes
~~~~~~~~~~~~~

* Support passing :class:`~capytaine.bodies.bodies.FloatingBody` or :class:`~capytaine.post_pro.free_surfaces.FreeSurface` objects to post-processing methods such as :meth:`~capytaine.bem.solver.BEMSolver.compute_potential` and :meth:`~capytaine.bem.solver.BEMSolver.compute_free_surface_elevation`. (:pull:`379`)

* Add ``top_light_intensity`` optional arguments to :meth:`~capytaine.ui.vtk.animation.Animation.run` and :meth:`~capytaine.ui.vtk.animation.Animation.save` to illuminate the scene from top. (:pull:`380`)

* Clean up ``__str__`` and ``__repr__`` representation of many objects. Also ``rich.print`` now return even nicer representations. (:pull:`384`)

* Always automatically compute and store the ``excitation_force`` next to the ``Froude_Krylov_force`` and ``diffraction_force`` in the dataset (:pull:`406`).

* Computing the RAO with :func:`~capytaine.post_pro.rao.rao` is not restricted to a single wave direction (or a single value of any other extra parameter) at the time anymore. (:issue:`405` and :pull:`406`)

* New computation of quadrature schemes without relying on Quadpy. (:pull:`416`)

* Add a new function :func:`~capytaine.io.legacy.run_cal_file` to solve the problems defined by a Nemoh.cal file, exactly as the command-line interface is doing (:pull:`422`).

* Rephrase mesh resolution warnings and group several of them together in a single warning. (:pull:`423`)

* Add block-Jacobi/coarse-correction preconditioner for large arrays of bodies. (:pull:`436`)

* The tabulation is saved on disk in a cache directory instead of being recomputed at each initialization of the solver. (:pull:`454`)

* Add a ``faces_max_radius`` argument to the predefined geometries from :mod:`~capytaine.meshes.predefined` to set up the resolution by giving a length scale for the panels (:pull:`459`).

* Automatically clip the mesh (and display a warning) when a problem is initialized with a mesh above the free surface or below the sea bottom (:pull:`486`).

Bug fixes
~~~~~~~~~

* When initializing a body with a mesh having degenerate panels, the initialization of the dofs used to happen before the degenerate panels were removed, leading to an inconsistency between the number of panels in the mesh and in the dof definition. (:issue:`367` and :pull:`375`)

* Fix the single precision Green function (:code:`cpt.Delhommeau(floating_point_precision="float32")`) that was broken in v2.0. (:issue:`377` and :pull:`378`)

* Update the BEMIO import feature to work with Pandas 2.0 and output periods as now done in Capytaine 2.0. A version of BEMIO that works in recent version of Python and Numpy can be found at https://github.com/mancellin/bemio. (:pull:`381`)

* Fix :meth:`~capytaine.bem.solver.BEMSolver.compute_pressure` that was broken. (:pull:`394`)

* Fix error message when computing hydrostatic stiffness of non-neutrally-buoyant body that is not a single rigid body. (:issue:`413` and :pull:`414`)

* Fix bug causing the quadrature method of a mesh to be forgotten when the mesh was put in a body. ``quadrature_method`` can now be passed as argument when initializing a new mesh. (:pull:`417`)

* The function :func:`~capytaine.io.mesh_loaders.load_mesh` more robustly detects filetype using file extension even when the file extension is not lowercase. (:pull:`441`)

* Fix bug with bodies translation or rotation when the rotation center or the center of mass had been defined as list or tuples instead of array (:pull:`472`).

Internals
~~~~~~~~~

* Add tentative build file for the Guix package manager (:pull:`339`).

* Fix badly named variables ``VSP2_SYM`` and ``VSP2_ANTISYM`` in libDelhommeau (:pull:`391`)

* Remove dependency to ``hypothesis`` for testing (:pull:`391`).

* Change how forces are stored in result objects. Added mass and radiation damping can now be queried with ``added_mass`` and ``radiation_damping`` and not only the plural forms that were used nowhere else in the code. (:pull:`393`)

* Use `nox <https://nox.thea.codes>`_ to test the code in isolated virtual environments. (:pull:`401`)

* Fortran source files are not included in wheel anymore (:pull:`360`).

* The ``delete_first_lru_cache`` decorator has been renamed :func:`~capytaine.tools.lru_cache.lru_cache_with_strict_maxsize` and now supports keyword arguments in the memoized function (:pull:`442`).

* Fix Xarray future warning about `Dataset.dims` (:issue:`450` and :pull:`451`).

* Improve some warnings and error messages.

-------------------------------
New in version 2.0 (2023-06-21)
-------------------------------

Major changes
~~~~~~~~~~~~~

* User can specify a period, a wavelength or a wavenumber instead of an angular frequency :code:`omega` when setting up a problem or a test matrix. If several types of frequency data are provided, an error is raised (:pull:`283`).

* **Breaking** The normalization of radiation problems has been changed to use the same normalization as diffraction problems. Added mass and radiation dampings are unchanged, but other outputs of radiation problem (free surface elevation, kochin functions, etc.) may differ from previous version by a factor :math:`-j \omega`. (:issue:`173` and :pull:`348`)

* **Breaking** The above two points interfered with the handling of :math:`\omega = 0` and :math:`\omega = \infty` cases. They have been temporarily disabled and will return in a future release.

* Add methods :meth:`~capytaine.bem.solver.BEMSolver.compute_potential`, :meth:`~capytaine.bem.solver.BEMSolver.compute_velocity` and :meth:`~capytaine.bem.solver.BEMSolver.compute_free_surface_elevation` and :meth:`~capytaine.bem.solver.BEMSolver.compute_pressure` to compute the value of some fields in the domain in post-processing. Their signature has been uniformized with the :func:`~capytaine.bem.airy_waves.airy_waves_potential` and :func:`~capytaine.bem.airy_waves.airy_waves_velocity` functions (:pull:`288`, :pull:`326`). New functions :func:`~capytaine.bem.airy_waves.airy_waves_free_surface_elevation` and :func:`~capytaine.bem.airy_waves.airy_waves_pressure` have also been added (:pull:`293`).

* **Breaking** The problems can now be initialized by setting a ``water_depth`` instead of the ``sea_bottom`` (which is still available for user-facing functions). This change is meant to uniformize notations in the code and use ``water_depth`` wherever possible (:pull:`340`). Besides the ``sea_bottom`` argument of many internal routines has been completely replaced by ``water_depth``. Migrating then requires changing the sign of the float (:pull:`347`).

* Add Github Actions workflow to build wheels. Precompiled packages will now be available with ``pip`` and not only with ``conda``.

Minor changes
~~~~~~~~~~~~~

* Support the new format of `Nemoh.cal` file from Nemoh v3 (:issue:`278` and :pull:`280`).

* **Breaking** Remove the :code:`convention` parameter to compute excitation force with WAMIT conventions (:issue:`133` and :pull:`281`).
  Changing the convention to compare the outputs of different codes is better done by a dedicated software such as `BEMRosetta <https://github.com/BEMRosetta/BEMRosetta>`_ or `BEMIO <https://wec-sim.github.io/bemio/>`_.

* Add nicer display for Capytaine objects in IPython shell (:issue:`227` and :pull:`287`).

* Support exporting hydrostatics data in original Nemoh-format files - :code:`Hydrostatics.dat` and :code:`KH.dat` (:pull:`285`).

* Add nicer display for Capytaine objects in IPython shell (:issue:`227` and :pull:`287`)

* Add functions :func:`~capytaine.io.mesh_loaders.load_PNL` and :func:`~capytaine.io.mesh_writers.write_PNL` to load and write meshes in HAMS ``.pnl`` format (:pull:`289`).

* **Breaking** Remove ``cpt.Nemoh()`` class that was replaced by :class:`~capytaine.bem.solver.BEMSolver` in version 1.1 (:pull:`291`)

* **Breaking** Remove ``full_body`` attribute from :class:`~capytaine.bodies.bodies.FloatingBody` that used to keep a copy of the body before clipping in-place (:pull:`302`).

* **Breaking** Remove ``dimensionless_wavenumber`` and ``dimensionless_omega`` attributes from :class:`~capytaine.bem.problems_and_results.LinearPotentialFlowProblem` as they are not used in the code and can be easily recomputed by users if necessary (:pull:`306`).

* Add :meth:`~capytaine.bodies.bodies.FloatingBody.minimal_computable_wavelength` to estimate the wavelengths computable with the mesh resolution (:pull:`341`).

* Slightly increase default tabulation size to avoid some high-frequency issues such as :issue:`157` (:pull:`353`).

Bug fixes
~~~~~~~~~

* Fix :meth:`~capytaine.meshes.collections.CollectionOfMeshes.immersed_part` (:pull:`307`).

* :meth:`~capytaine.bodies.bodies.FloatingBody.compute_hydrostatics` used to fail for non-rigid bodies because it could not compute the rigid-body inertia.
  The rigid-body inertia is now just skipped for bodies with no rigid-body dofs (:pull:`308`).

* Reduce the default tolerance of the mesh clipper for points just above the free surface (:issue:`320` and :pull:`322`).

* Convert ``center_of_mass`` and ``rotation_center`` to arrays in :class:`~capytaine.bodies.bodies.FloatingBody` constructor to avoid a few issues (:issue:`319` and :pull:`325`).

* Fix bug (leading to either ``RuntimeError`` or wrong output) when clipping with plane that does not contain the origin. (:pull:`344`)

* Instances of :class:`~capytaine.bem.solver.BEMSolver` initialized with default parameters do not share the same engine, hence they do not share the same cache. This minor issue was causing minor interferences in some benchmarks (:issue:`295` and :pull:`350`).

Internals
~~~~~~~~~

* Major update of the compilation toolchain because of the upcoming deprecation of ``numpy.distutils``. Capytaine is now built with ``meson-python``.

* The method :meth:`~capytaine.green_functions.delhommeau.Delhommeau.evaluate` (and its counterparts for other Green functions) now accepts a list of points as first argument instead of a mesh. It has now an optional boolean argument ``early_dot_product`` to return the integrals of the gradient of the Green function and not only the normal derivative (:pull:`288`).

* Remove warnings due to 0/0 divisions in :func:`~capytaine.meshes.properties.compute_faces_properties` (:pull:`310`)

* **Breaking** Remove unused and undocumented code about meshes, including ``mesh.min_edge_length``, ``mesh.mean_edge_length``, ``mesh.max_edge_length``, ``mesh.get_surface_integrals``, ``mesh.volume``, ``mesh.vv``, ``mesh.vf``, ``mesh.ff``, ``mesh.boundaries``, ``mesh.nb_boundaries``, ``compute_faces_integrals``, ``SingleFace``. (:pull:`334`)

* Add analytics to the documentation using `Plausible.io <https://plausible.io>`_ (:pull:`290`).

-------------------------------
New in version 1.5 (2022-12-13)
-------------------------------

Major changes
~~~~~~~~~~~~~

* The :class:`~capytaine.green_functions.delhommeau.XieDelhommeau` implementation of the Green function has been improved.
  The implementation used to be almost the same as the default :class:`~capytaine.green_functions.delhommeau.Delhommeau` method.
  A missing key element has been added and the :class:`~capytaine.green_functions.delhommeau.XieDelhommeau` is now actually more accurate near the free surface.
  (:pull:`180` and :pull:`216`)

* New default linear solver :class:`~capytaine.matrices.linear_solvers.LUSolverWithCache`: the LU decomposition of the matrix is now cached to be reused for other similar problems, diminishing the total computation time up to 40%. (:pull:`235`)

* New functions to generate simple geometric meshes have been implemented in :code:`capytaine.meshes.predefined`. They are similar to the former geometric bodies (:class:`~capytaine.bodies.predefined.sphere.Sphere`, :class:`~capytaine.bodies.predefined.sphere.HorizontalCylinder`, etc.), except that they return a mesh and do not create a :code:`FloatingBody`. The geometric body classes are considered deprecated, although they should still work as expected. (:pull:`233`)

* Changed the behavior of :meth:`~capytaine.bodies.bodies.FloatingBody.compute_hydrostatics`. The mesh is not silently modified anymore. The stiffness and inertia matrices are stored in the body for inclusion in the output dataset. The inertia matrix is now computed on the full mesh (:issue:`197`, :issue:`249`, :issue:`258` and :pull:`262`).

Minor changes
~~~~~~~~~~~~~

* Add :code:`floating_point_precision` argument to :meth:`~capytaine.green_functions.delhommeau.Delhommeau` and :meth:`~capytaine.green_functions.delhommeau.XieDelhommeau` that accepts either :code:`"float32"` for single precision computations or :code:`"float64"` for double precision computations (the latter is the default). (:pull:`224`).

* Passing the argument :code:`tabulation_nr=0` or :code:`tabulation_nz=0` to :class:`~capytaine.green_functions.delhommeau.Delhommeau`
  or :class:`~capytaine.green_functions.delhommeau.XieDelhommeau` now allows to run the code without interpolating the Green function
  from a precomputed tabulation. This is meant as a tools for benchmarks and validation, since it decreases the performance of the code
  for often no accuracy gain. (:pull:`229`)

* :func:`~capytaine.io.mesh_loaders.load_mesh` is now exported by the main namespace: :code:`from capytaine import load_mesh`.
  The documentation has been changed to recommend the use of this function instead of :meth:`~capytaine.bodies.bodies.FloatingBody.from_file`.
  (:pull:`231`)

* When initializing a :code:`FloatingBody`, one can now pass directly a mesh object from :code:`meshio`.
  The documentation has been changed to recommend this approach instead of :meth:`~capytaine.bodies.bodies.FloatingBody.from_meshio`.
  (:issue:`259` and :pull:`261`)

* When joining two bodies as e.g. :code:`body1 + body2`, some hydrostatic properties are passed to the resulting body:
  if all the bodies have hydrostatic stiffness matrices or inertia matrices defined,
  then they are assigned to the joined body as a larger block diagonal matrix (:pull:`243`).

* Add :meth:`~capytaine.bodies.bodies.FloatingBody.immersed_part` method to clip the body without modifying it in place (:pull:`244`).

* Add :func:`~capytaine.rigid_body_dofs` method returning a placeholder that can be given at the creation of :class:`~capytaine.bodies.bodies.FloatingBody` to initialize the six rigid body dofs (:pull:`245`).

* Custom classes from the :code:`capytaine.matrices` module storing block matrices or data-sparse matrices
  can be transformed into full Numpy arrays with :code:`np.array(...)` (:pull:`99`)

* Add :code:`Dockerfile` and instructions to install with Docker (:pull:`137`)

* Add optional arguments to :func:`~capytaine.io.meshes_writers.write_GDF` to write parameters :code:`ulen, grav, isx, isy` to the mesh file (:pull:`241`)

* Fix bug with MED mesh file loading (:issue:`247` and :pull:`250`).

* Several surface integrals properties of :code:`FloatingBodies` are also defined on meshes, such as :code:`volume` or :code:`center_of_buoyancy` (pull:`263`).

Internals
~~~~~~~~~

* The integration of the pressure on the mesh of the body was implemented twice independently. It has been factored out in :meth:`~capytaine.bodies.bodies.FloatingBody.integrate_pressure` (:pull:`218`)

* `__rmatmul__` has been implemented for low rank matrices (:pull:`222`).

* New implementation of the GDF mesh file reader :func:`~capytaine.io.meshes_loaders.load_GDF` (:pull:`241`)

---------------------------------
New in version 1.4.2 (2022-10-03)
---------------------------------

Bug fixes
~~~~~~~~~

* Raise error message when calling :meth:`~capytaine.bodies.bodies.FloatingBody.compute_hydrostatics()` without a center of mass defined (:pull:`207`).

* Fix bug when cropping body with a dof defined manually as a list of tuples (:issue:`204` and :pull:`206`).

Documentation
~~~~~~~~~~~~~

* Miscellaneous improvements of the documentation (:pull:`205`, :pull:`211`, :pull:`219`)

* Clean up and fix animation example in the cookbook (:pull:`213`).

* The warning message for insufficient mesh resolution appears earlier and has been reworded to be clearer (:pull:`217`).

Internals
~~~~~~~~~

* Replace the Fortran core by a git submodule pointing to `libDelhommeau <https://github.com/capytaine/libDelhommeau/>`_ (:pull:`208`).
  Future developments of the Green function will take place there.

* Move from Travis CI to Github Actions for continuous integration (:pull:`209`)

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
  Nothing changes when users explicitly choose a linear solver. (:pull:`171`)

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
