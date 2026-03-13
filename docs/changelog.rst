=========
Changelog
=========

.. contents::
   :local:
   :depth: 1
   :backlinks: none

Older changelogs can be found on this page: :doc:`older_changelog`.

-------------------------------
New in version 3.0 (2026-??-??)
-------------------------------

Major changes
~~~~~~~~~~~~~

* Add :class:`~capytaine.bodies.multibodies.Multibody` meant to represent a
  multibody system. For hydrodynamics, this is equivalent to the previous
  behavior of coalescing bodies together. For hydrostatics, the new class
  is slightly more powerful, for instance by being able to keep track of
  several center of buoyancy and center of mass. (:pull:`822`)

  Joining bodies with :meth:`~capytaine.bodies.bodies.FloatingBody.join_bodies`
  or ``+`` now creates a :class:`~capytaine.bodies.multibodies.Multibody` instance. It can be converted back to a
  :class:`~capytaine.bodies.bodies.FloatingBody` instance with::

    both = body_1 + body_2  # `both` is now a Multibody
    both = (body_1 + body_2).as_FloatingBody()  # Recover former behavior of joining FloatingBody with a FloatingBody

* New internal data model for rigid body dofs with the classes
  :class:`~capytaine.bodies.dofs.TranslationDof` and
  :class:`~capytaine.bodies.dofs.RotationDof`. (:pull:`838`)
  This should ensure that rigid dofs are detected and treated as such, without
  relying only on their name.
  This is relevant for multiple bodies and articulated bodies when computing:
    * hydrostatics, where the exact hydrostatic stiffness formula for rigid body dofs can be used instead of the approximation for generalized dofs.
    * forward speed, where the m-term is currently only implemented for rigid dofs.


Minor changes
~~~~~~~~~~~~~

* Add option :code:`'lu_decomposition_with_overwrite'` for the :code:`linear_solver` of :class:`~capytaine.bem.engines.BasicMatrixEngine`, which reduces the RAM usage of the solver (:pull:`775`).

* Velocity in the fluid can be post-processed in the limit frequencies (:math:`\omega = 0` or :math:`\omega = \infty`). Divide it by :math:`\omega` to have a finite value (:pull:`777`).

* Display in the log the RAM usage estimation before a batch resolution and the measured RAM usage at the end of the resolution (:pull:`784`)

* **Breaking** When building the dataset in :meth:`~cpt.bem.solver.fill_dataset`, the previous ``kochin`` attribute has been renamed to ``kochin_radiation``
  to be consistent with the existing ``kochin_diffraction`` attribute.

* **Breaking** Remove the geometric body classes ``cpt.Sphere()``,
  ``cpt.VerticalCylinder()``, ``cpt.HorizontalCylinder()``, ``cpt.Disk()``,
  ``cpt.Rectangle()``, ``cpt.RectangularParallelepiped()``, that were marked as
  deprecated since version 2.0. Consider using instead::

    cpt.FloatingBody(mesh=cpt.mesh_sphere(...), ...)

  or something equivalent, separating the geometric mesh generation from the
  floating body definition.

* :func:`~capytaine.bodies.dofs.rigid_body_dofs` now instantiate the new dof class instead of a placeholder.
  It now has an additional input argument ``only`` to have only some of the six rigid body dofs.
  Also a rotation center should be passed, because the properties of the body
  (e.g. `center_of_mass`) cannot be accessed at that stage. (:pull:`838`)

* Add function :func:`~capytaine.post_pro.mean_drift_force.far_field_mean_drift_force` to compute the horizontal mean drift forces using far field formulation.
  Only the single-direction second-order term is currently implemented.

Bug fixes
~~~~~~~~~

* Fix type of the right-hand-side of the linear solver when using option :code:`floating_point_precision = 'float32'` in :class:`~capytaine.green_functions.delhommeau.Delhommeau`.
  As a consequence, the whole computation is done in single precision and the RAM usage is lower as expected. (:pull:`774`)

* Fix the timer for parallel resolution, i.e. when :code:`n_jobs` is greater than 1. Now the durations for each process are displayed. (:pull:`782`)

* Hydrostatics methods better take into account the free surface even when
  the mesh has not been clipped yet.
  Internally, the mesh is clipped before computing methods such as
  :meth:`~cpt.bodies.bodies.FloatingBody.disp_volume` or
  :meth:`~cpt.bodies.bodies.FloatingBody.center_of_buoyancy`, such that the
  computed displaced volume is always actually the volume below $z=0$ and the
  center of buoyancy is always below the free surface.
  The inertia matrix always uses the displaced water mass for the mass and
  compute the inertia moments on the full shape if the full mesh is provided.
  (:pull:`794`)

* Fix the frequency type in the dimensions of the dataset returned by :meth:`~cpt.io.xarray.kochin_data_array`. Previously, the dimension was always named ``omega``;
  it is now named ``omega``, ``freq``, ``period``,  ``wavenumber`` or ``wavelength`` depending on the user settings.

* Fix a bug when computing the Kochin function on a mesh with a lid. The Kochin function now also takes the faces on the lid into account. (:issue:`833`)

* Fix WAMIT output in ``.1`` file where 0 and infinite frequency added mass where interverted. (:pull:`855`)

Internals
~~~~~~~~~

* **Breaking** The ``green_function`` is not an attribute of the :class:`~capytaine.bem.solver.BEMSolver` anymore, but of the engine.
  The motivation is that not all engines can be made compatible with all Green function implementations (although the builtins one are).
  The possibility to call ``BEMSolver(green_function=...)`` is kept as a convenient shortcut to ``BEMSolver(engine=BasicMatrixEngine(green_function=...))``.
  Calls to ``BEMSolver(green_function=..., engine=...)`` now raise an error. (:pull:`752`)
  Post-processing new requires the implementation of the methods ``build_S_matrix`` and ``build_fullK_matrix`` by the engine (:pull:`753`)

* New implementation of the block symmetric matrices for mesh symmetry, now
  used by :class:`~capytaine.bem.engines.BasicMatrixEngine` (:pull:`754`).

* Rafactor of the :class:`~capytaine.bem.engines.BasicMatrixEngine` to make the
  caching more straightforward and improve its interaction with LU
  decomposition and symmetries. (:pull:`755`)

* The whole ``matrices`` module as well as the corresponding engines in
  ``bem/engines.py`` and ``tool/lru_caches.py`` have been removed from
  Capytaine. For compatibility, they will remain accessible from a separate
  package. (:pull:`757`, :pull:`765`)

* Instead of creating a ``CollectionOfMeshes``, now ``join_bodies`` merges the
  meshes of the bodies together. The legacy behavior is still available from
  Fakeblocks ``join_bodies``. (:pull:`779`)

* The new ``Mesh`` class has a ``faces_metadata`` attribute storing fields
  defined on each faces of the mesh. When faces are added or removed, the
  metadata are automatically updated accordingly. (:pull:`791`)

* Move hydrostatics routines in a dedicated module and rewrite corresponding tests (:pull:`794`)

* Refactor the implementation of the timer to make it easier to include more steps (:pull:`809`)

* Parameters ``free_surface``, ``water_depth`` and ``wavenumber`` are always
  keyword arguments in the engine and the Green function (:pull:`812`)

* **Breaking** :meth:`~capytaine.bodies.FloatingBody.add_rotation_dof` takes
  arguments called ``rotation_center`` and ``direction`` instead of an ``Axis``
  object. Also the geometric center is not used anymore as a fallback value for
  ``rotation_center``.

* The S matrix can now be computed even if the K matrix is not defined,
  for example when evaluating the :meth:`~capytaine.bem.solver.compute_free_surface_elevation` along the waterline.

---------------------------------
New in version 2.3.1 (2025-10-14)
---------------------------------

Bug fix
~~~~~~~

* Fix **major bug of version 2.3** where the resolution of problem with **both mesh
  symmetries and a lid** for irregular frequencies removal returned wrong values.
  (:issue:`761`)

* Fix issue where in-place transformation of a ``FloatingBody`` (such as
  ``body.keep_immersed_part()`` or ``body.translate(...)``) were sometimes not
  taken into account. In-place transformation are not recommended and might be
  removed in a future version, use the versions returning new objects as seen
  in the documentation (e.g. ``body.immersed_part()`` and
  ``body.translated(...)``).

* If loading the tabulation from the file fails, then the tabulation is
  recomputed (`Issue 739
  <https://github.com/capytaine/capytaine/issues/739#issuecomment-3190735343>`_)

Internals
~~~~~~~~~

- The source code moved from ``capytaine`` to ``src/capytaine`` in the main
  repository to avoid importing the local folder instead of the installed
  version (:issue:`395` and :pull:`749`).

- Replace development dependencies in ``editable_install_requirements.txt`` and
  ``[project.optional-dependencies]`` with ``[dependency-groups]``
  (:pull:`750`).


-------------------------------
New in version 2.3 (2025-07-17)
-------------------------------

Major change
~~~~~~~~~~~~

* The implementations of the **Green function from HAMS** are now included in Capytaine:

  * The infinite depth version from [Liang, Wu, Noblesse, 2018] is :class:`~capytaine.green_functions.hams.LiangWuNoblesseGF` (:pull:`617`),
  * The finite depth version from [Liu et al., 2018] is :class:`~capytaine.green_functions.hams.FinGreen3D` (:pull:`647`),
  * The class :class:`~capytaine.green_functions.hams.HAMS_GF` is a thin wrapper using one or the other method above depending of the water depth (:pull:`658`).

  They can be passed to Capytaine's solver as follows::

    solver = cpt.BEMSolver(green_function=cpt.HAMS_GF())

  Please cite the corresponding papers if you use them in a scientific publication (see the :doc:`citing` page).

* Revamp of default **finite depth Green function** implementation.

  * The new implementation should better handle panels on or near the free surface and have the right asymptotic consistency with the infinite depth method when depth goes to infinity.
    The legacy behavior of previous versions is still available by setting the parameter :code:`finite_depth_method` added to :class:`~capytaine.green_functions.delhommeau.Delhommeau` to :code:`finite_depth_method="legacy"`, while the better behavior is used by default. (:pull:`654` and :pull:`656`)
  * The Prony decomposition is now done in Python and its failure (typically for :math:`kh < 0.1`) raises an error instead of returning wrong values.
    This behavior is controlled by the :code:`finite_depth_prony_decomposition_method` parameter of :class:`~capytaine.green_functions.delhommeau.Delhommeau`, which is now :code:`"python"` by default. (:pull:`675`)
  * Infinite frequency is now supported in finite depth (zero frequency is still not and returns the same error as other finite depth low-frequency cases). (:pull:`703`)

* Do not interrupt a batch of resolutions when one of them fails.
  Instead the exception is displayed in the log and the results are replaced by a :class:`~capytaine.bem.problems_and_results.FailedDiffractionResult` or :class:`~capytaine.bem.problems_and_results.FailedRadiationResult`. The output dataset is filled with a ``NaN`` value for these parameters. (:pull:`678`)
  Diffraction problems with zero or infinite frequencies used to have a special treatment to be run with a batch resolution despite raising an error when run alone, they have been reworked to use the same design as other failing resolutions. (:pull:`719`)

* The Boundary Integral Equation (``method`` keyword argument) used to solve the problem can now be specified when initializing a solver and will then be use for all resolution with this solver. This general setting can be over overridden by using the ``method`` argument when solving::

    solver = cpt.BEMSolver(method="direct")  # That is new and recommended
    solver.solve(problem, method="direct")  # That is still possible and override the above setting.

The method is also saved in the metadata of the results with the other parameters of the solver (whether it was defined when initializing the solver or later). (:pull:`686`)

* Add :func:`~capytaine.io.xarray.export_dataset` method to more conveniently export a dataset to NetCDF or other formats (:pull:`690`).

Minor change
~~~~~~~~~~~~

* Add optional :code:`freq` argument (frequency in Hz) for problem set up and output.

* Add :func:`~capytaine.io.xarray.assemble_dataframe` which collect results into a Pandas DataFrame (this was already done internally in :func:`~capytaine.io.xarray.assemble_dataset`) (:pull:`677`).
  Also add :func:`~capytaine.io.xarray.assemble_matrices` function which is a simplified version of :func:`~capytaine.io.xarray.assemble_dataset` without metadata, meant to be used mostly for teaching. (:pull:`643`)

* The environment variable ``CAPYTAINE_PROGRESS_BAR`` can be used to disable globally the display of a progress bar when solving problems. This is meant mostly for testing environments and CI. (:pull:`646`)

* Add ``timer`` attribute to :class:`~capytaine.bem.solver.BEMSolver` storing the time spent in each steps of the resolution. Summary can be accessed by :meth:`~capytaine.bem.solver.BEMSolver.timer_summary`. (:pull:`674`)

* Add :func:`~capytaine.io.wamit.export_to_wamit` as a unified interface to export hydrodynamic results to WAMIT-compatible files. (:pull:`714`)


Bug fixes
~~~~~~~~~

* Properly use ``progress_bar`` argument in :func:`~capytaine.bem.solver.fill_dataset` to disable progress bar.

* Always remove degenerate faces after clipping (:issue:`620` and :pull:`624`).

* Fix missing geometric center in legacy predefined body :class:`~capytaine.bodies.predefined.rectangles.ReflectionSymmetricMesh`. It was causing inconsistent definition of dofs with respect to earlier versions. (:pull:`625`)

* Fix Python implementation of the Prony decomposition for the finite depth Green function, which is now the default. (:pull:`621`). Move some code of its code to the :mod:`~capytaine.tools.prony_decomposition` module. (:pull:`649`)

* After joining several bodies, editing the mesh of one of the components does not affect the joined body anymore (:issue:`660` and :pull:`662`:).

* Check the consistency of the dofs with the mesh and raises ``ValueError`` when an inconsistency is detected (:pull:`663`).

* Fix error when removing all the faces from a symmetric mesh (:pull:`668`)

* Add safeguard if a custom linear solver returns a result vector of wrong shape (e.g. column instead of row) (:pull:`670`)

* Fix loading BEMIO datasets from Nemoh (:pull:`681`)

* Fix computing zero and infinite frequency radiation problems with a lid for irregular frequencies removal (:issue:`704` and :pull:`708`)

* Fix solving :class:`~capytaine.bem.problems_and_results.LinearPotentialFlowProblem` directly.

* Fix missing variable attributes for main frequency variable (:issue:`702` and :pull:`717`)

* Trying to generate a lid over a purely vertical mesh does not raise an error anymore (:issue:`625`).

* When the hull mesh and the lid mesh are both symmetric with the same reflection plane, the symmetry is not lost anymore when solving the BEM problem.
  Also ``generate_lid`` and ``extract_lid`` should now work with reflection symmetric meshes without losing the symmetry. (:issue:`527`, :pull:`667`, :pull:`720`).


Internals
~~~~~~~~~

* Major refactoring of the Fortran core, including its interface in Python:

  * Add ``interface.f90`` Fortran file to group some routines used only for wrapping the Fortran core. (:pull:`612`)

  * Add :meth:`~capytaine.green_functions.delhommeau.Delhommeau.all_tabulation_parameters` to make it easier to test Fortran core from Python (:pull:`648`)

  * Refactor implementation of Delhommeau's finite depth Green function to compute all the frequency-independant Rankine terms at the same time (for future caching) (:pull:`652`)

  * The main interface to the Fortran core ``build_matrices`` does not take ``coeffs`` and ``same_body`` inputs anymore.
    The role of the former is played by ``gf_singularities`` and ``wavenumber``.
    The diagonal term added by the latter is now added independently.
    (:pull:`701`)

* NaN values are not striped out of output data (:pull:`676`)

* Define a :class:`~capytaine.meshes.mesh_like_protocol.MeshLike` protocol that classes implementing a mesh should follow. Also ensure that :class:`~capytaine.meshes.meshes.Mesh` and :class:`~capytaine.meshes.collections.CollectionOfMeshes` follow it. (:pull:`667`)

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
  Meshes files can list points in either 3×4 × ``nb_faces`` or a 12 × ``nb_faces`` format.
  (:issue:`540` and :pull:`585`).

* When filling a test matrix with both diffraction problems and radiation
  problems, zero and infinite frequencies can now be provided. (Previously, the
  computation was failing because these frequencies are not defined for
  diffraction problems.) (:pull:`587`)

* Radiation damping at infinite frequency is now zero instead of infinity.
  When forward speed is non-zero, added mass and radiation dampings at zero encounter frequency are NaN.
  (:pull:`588`)

* User does not need to import ``pyplot`` themself before running :meth:`~capytaine.meshes.meshes.Mesh.show_matplotlib()` (:pull:`592`)

* Fixes usage of ``ReflectionSymmetricMesh`` with direct solver (:issue:`593` and :pull:`594`).

* Do not recompute the same
  :meth:`~capytaine.bodies.bodies.FloatingBody.first_irregular_frequency_estimate`
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

* Meshio objects can be directly passed to :func:`~capytaine.io.mesh_loaders.load_mesh` to get a Capytaine mesh (:pull:`555`).

* Load gmsh v4 format .msh file using :func:`~capytaine.io.mesh_loaders.load_mesh` (when meshio is installed) (:pull:`556`)


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
