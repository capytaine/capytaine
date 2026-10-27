=========
Changelog
=========

This page is for changes introduced in versions 3.x of Capytaine.
For older changelogs, see:

.. toctree::
   :maxdepth: 1

   changelog_from_2_to_3
   changelog_before_2

-------------------------------
New in version 3.1 (2026-??-??)
-------------------------------

Major changes
~~~~~~~~~~~~~

* Add a :class:`~capytaine.bodies.dofs.CustomDof` to define any generalized body mode.
  The legacy style of defining a dof with an array of motions allowed that already, but this new interface also allows optionally to give the gradient (actually jacobian) of the motion which is useful for hydrostatics, forward speed and second-order forces.
  Passing the ``divergence`` parameter to :meth:`~capytaine.bodies.hydrostatics.compute_hydrostatic_stiffness` is discouraged in favor of this new construct.
  (:pull:`939`)

* :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` has a new keyword argument ``keep_details`` (default: ``False``).
  When set to ``True``, the pressure field on the hull is kept for each problem and automatically included in the returned dataset as the data variables ``pressure_diffraction``, ``pressure_Froude_Krylov`` and ``pressure_radiation``.
  The source distribution of the indirect method is also stored, as ``sources_diffraction`` and ``sources_radiation`` (on the hull, with dimension ``hull_face``) and ``sources_lid_diffraction`` and ``sources_lid_radiation`` (on the lid, if any, with dimension ``lid_face``).
  With the direct method, the potential is stored instead, as ``potential_diffraction``, ``potential_radiation``, ``potential_lid_diffraction`` and ``potential_lid_radiation``.
  Results of the direct and indirect methods solved with ``keep_details=True`` cannot be assembled in the same dataset (they can with ``keep_details=False``).
  See :doc:`user_manual/post_pro` for more details.

* The ``mesh=True`` option of :func:`~capytaine.io.xarray.assemble_dataset` and :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` now stores the full mesh of the hull (and lid, if any) and the dofs of the body in the dataset, as the data variables ``mesh_vertices``, ``mesh_faces_center``, ``lid_mesh_vertices``, ``lid_mesh_faces_center``, ``dof_motions`` and ``dof_gradient_of_motions``, instead of just the number of faces.
  See :doc:`user_manual/post_pro` for more details.

* Add :func:`~capytaine.post_pro.mean_drift_force.near_field_mean_drift_force` to compute the mean drift force, for any dofs of a :class:`~capytaine.bodies.bodies.FloatingBody` or a :class:`~capytaine.bodies.multibodies.Multibody` (rigid body motions or generalized modes with a gradient of motion). (:pull:`871`)
  It is called as ``near_field_mean_drift_force(first_order_motion_amplitude, dataset, body)``, where ``dataset`` has been computed with the indirect method and ``keep_details=True`` (e.g. ``solver.fill_dataset(test_matrix, body, keep_details=True)``): the first-order fields are evaluated from the source distributions stored in the dataset, without solving the BEM problems again.
  The force is computed for the frequencies and wave directions common to ``first_order_motion_amplitude`` and ``dataset``, and for the dofs common to them and ``body`` (a warning is logged when some are left out).
  This function returns a dataset with a ``near_field_mean_drift_force`` variable of shape ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_dofs)``.
  Optional parameter ``output_pressure`` allows to return the mean drift pressure field on the hull as ``pressure_hull_mean_drift``, of shape ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_hull_faces)``,
  as well as the first-order fields on the waterline from which the waterline term of the force is computed: ``pressure_waterline_total`` and ``relative_elevation_waterline``, of shape ``(nb_freq, nb_wave_direction, nb_waterline_segments)``, and the waterline segments ``mesh_waterline``.
  The near field mean drift force can be exported to WAMIT's `.9` file format. (:pull:`909`)

Bug fixes
~~~~~~~~~

* For a mesh without waterline (e.g. a fully submerged body), ``edges_waterline`` and ``faces_waterline`` are now empty 2D arrays (of shape ``(0, 2)`` and ``(0, 4)``) instead of 1D arrays, so that the computations on the waterline do not fail.

* Fix the conversion of a :class:`~capytaine.tools.block_circulant_matrices.BlockDiagonalMatrix` to a NumPy array, that was raising an ``AttributeError``. (:pull:`958`)

* Fix computation of hydrostatic stiffness on multibodies treated as generalized dofs as in previous version. (:pull:`936`)

* Fix double layer term when computing the value at field points without a free surface.

* The GDF reader of :func:`~capytaine.meshes.io.load_mesh` now follows the free format of WAMIT: the coordinates of a panel can be spread over any number of lines, extra tokens after them (e.g. a panel number) and lines after the ``NPAN`` panels (e.g. ``NPAND``) are ignored, and commas and ``D`` exponents are accepted.
  A warning is emitted for the dipole panels of the WAMIT Version 6 format, which are not supported, and an error is raised for walls (negative ``ISX`` or ``ISY``) and higher-order geometries (``IGDEF`` other than 0).

* Fix a memory leak when solving with a mesh having two nested symmetries (two reflection planes or a dihedral symmetry): the influence matrices of all the previously solved frequencies were kept in memory.

* Fix the hash of Green functions, which depended on the identity of the object instead of on its settings: two identical Green functions now have the same hash.

* Fix :meth:`~capytaine.meshes.meshes.Mesh.extract_faces` on a mesh with faces metadata, that was failing. It now also accepts a slice or a boolean mask to select the faces.

Internals
~~~~~~~~~

* Methods about LU decomposition of numpy arrays have moved from :mod:`~capytaine.tools.block_circulant_matrices` to the new module :mod:`~capytaine.bem.linear_solvers`.
  :func:`~capytaine.bem.linear_solvers.lu_decompose` is now a :func:`functools.singledispatch` function,
  and the LU decompositions are subclasses of the new abstract class :class:`~capytaine.bem.linear_solvers.AbstractLUDecomposedMatrix`.
  The LU decomposition of a new type of matrices can be supported by the linear solver of the :class:`~capytaine.bem.engines.DefaultMatrixEngine` by registering it with ``lu_decompose.register``. (:pull:`954`)

* Add the dependency `array-api-compat <https://data-apis.org/array-api-compat/>`_ (and `array-api-strict <https://github.com/data-apis/array-api-strict>`_ to the test dependencies).
  The new module :mod:`capytaine.tools.array_backend` provides helpers (``is_array``, ``complex_dtype``, ``to_backend_of``, ``to_numpy``) to write code that works with any array library implementing the array API standard. (:pull:`956`)

* The :class:`~capytaine.bem.solver.BEMSolver` and the :class:`~capytaine.bem.engines.DefaultMatrixEngine` can now work with a Green function returning matrices from another array library than NumPy, as long as it implements the array API standard (using ``array-api-compat``) and an LU decomposition has been registered with :func:`~capytaine.bem.linear_solvers.lu_decompose`. (:pull:`957`)
  The results and the post-processing stay NumPy arrays.
  The subclasses of :class:`~capytaine.bem.linear_solvers.AbstractLUDecomposedMatrix` now have to define the method ``__array_namespace__`` and the attribute ``device``, as arrays of the array API standard do.
  The Green functions have new attributes ``matrices_array_backend`` and ``matrices_device`` to declare the array library and the device of their matrices (NumPy and ``"cpu"`` by default).
  The ``gmres`` linear solver raises ``NotImplementedError`` for matrices that are not NumPy arrays.
  Meshes with symmetries are supported: the block matrices of :mod:`capytaine.tools.block_circulant_matrices` are now implemented with the array API standard instead of NumPy. (:pull:`958`)
  They can be built from an array of shape ``(nb_blocks, n, m)`` or a list of blocks.
  The blocks of a :class:`~capytaine.tools.block_circulant_matrices.BlockCirculantMatrix` have to be arrays and not other block matrices anymore (use :class:`~capytaine.tools.block_circulant_matrices.NestedBlockCirculantMatrix` instead).

* The string representation (``__str__``, ``__repr__``) and the hash of Green functions are now implemented once in :class:`~capytaine.green_functions.abstract_green_function.AbstractGreenFunction`, driven by ``exportable_settings`` and a ``_default_parameters`` class attribute.

* **Rework the dof module**:

  * The rest of the code now only evaluates dofs on meshes.
    :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_motion` and
    :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_gradient_of_motion`
    take an optional argument ``location`` (``"faces_centers"`` by default,
    ``"quadrature_points"`` or ``"waterline"``).
    ``evaluate_motion_at_points`` and ``evaluate_gradient_of_motion_at_points`` are now considered
    the low-level definition of the dofs, to be implemented by most subclasses of
    :class:`~capytaine.bodies.dofs.AbstractDof`, except when reference to the mesh is required
    like for :class:`~capytaine.bodies.dofs.DofOnSubmesh`.
    User implemented subclasses of ``AbstractDof`` needs no change.

  * The implementations of :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_motion_at_points` and
    :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_gradient_of_motion_at_points`
    accept arrays of points of shape ``(..., 3)`` instead of only ``(nb_points, 3)``.

  * Remove the ``evaluate_motion_at_points`` and ``evaluate_gradient_of_motion_at_points``
    methods of :class:`~capytaine.bodies.dofs.DofOnSubmesh` which
    were indexing arbitrary points as if they were mesh faces. (:pull:`938`)

  * :class:`~capytaine.bodies.dofs.DofOnSubmesh` now only accepts a boolean array of size
    ``nb_faces`` to define its faces (not anymore a list of indices or a slice).

  * Add :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_second_order_motion` to the dofs,
    returning the second-order term of the expansion of the motion with respect to the amplitude of a pair of dofs.
    It is zero except for pairs of rotations around the same center, for which it follows the roll-pitch-yaw
    sequence of rotations of WAMIT (``R = R_roll R_pitch R_yaw``).
    Also add :meth:`~capytaine.bodies.dofs.AbstractDof.evaluate_gradient_of_second_order_motion`
    returning the gradient of the second-order motion.
    They take a location argument like the first order ones.
    For :class:`~capytaine.bodies.dofs.DofOnSubmesh`, it is zero outside of the faces common to both dofs
    (e.g. for two rotations of two different bodies of a multibody).

* Add a :code:`pressure_on_hull` attribute to
  :class:`~capytaine.bem.problems_and_results.LinearPotentialFlowResult`, equal
  to :code:`pressure` restricted to the faces of the hull mesh, excluding the
  lid used to remove irregular frequencies (if any).

* Speed up the test suite by about a factor of 2, by running it with
  ``OMP_NUM_THREADS=1`` (multithreading is mostly overhead on the small meshes of
  the tests) and by sharing the Green function instances between tests (instead
  of reloading the tabulation at each test). Helpers shared by several test
  files are gathered in ``pytest/capytaine_test_helpers.py``.

* Add :class:`~capytaine.tools.perturbation_expansion.PerturbationExpansion` to
  compute product of first and second order expansions, for the second-order
  forces. It supports batched fields (e.g. one matrix per face), accepts plain
  arrays and scalars as order-0 operands, and provides ``trace``, ``sum``,
  ``expand_dims``, ``squeeze``, ``det`` and ``inv``.
  Its subclass :class:`~capytaine.tools.perturbation_expansion.DifferenceFrequencyExpansion`
  holds complex amplitudes for several wave components at once, and computes the full
  difference-frequency QTF of an expression in one evaluation.
  First-order quantities are built with :meth:`~capytaine.tools.perturbation_expansion.DifferenceFrequencyExpansion.first_order`.
  :func:`~capytaine.post_pro.mean_drift_force.near_field_mean_drift_force` is now computed with :class:`~capytaine.tools.perturbation_expansion.DifferenceFrequencyExpansion`,
  from the unexpanded integral of the pressure on the moving body, for any dofs
  with a gradient of motion, with the help of
  :meth:`~capytaine.bodies.bodies.FloatingBody.integrate_pressure_on_moving_body`
  and
  :meth:`~capytaine.bodies.bodies.FloatingBody.integrate_pressure_on_waterline`.
  Add the properties ``middle_edges_waterline``, ``normals_edges_waterline`` and ``indices_faces_waterline`` to the meshes.

* Add the private methods ``BEMSolver._compute_potential_from_sources`` and ``BEMSolver._compute_potential_gradient_from_sources``, evaluating the potential or its gradient from bare source distributions (given with the mesh and the Green function parameters, without a result object).
  The sources can have leading batch axes, in which case the influence matrix is built only once and multiplied by all the source distributions at once.
  :meth:`~capytaine.bem.solver.BEMSolver.compute_potential`, :meth:`~capytaine.bem.solver.BEMSolver.compute_velocity`, :meth:`~capytaine.bem.solver.BEMSolver.compute_pressure` and :meth:`~capytaine.bem.solver.BEMSolver.compute_free_surface_elevation` are now thin wrappers around them.
  The ``LazyMatrix`` class now supports products with a 2D right-hand side.


-------------------------------
New in version 3.0 (2026-08-21)
-------------------------------

Version 3 includes a major rewrite of several internal modules, most notably
the ``meshes`` module which has been rewritten from scratch.

About the new mesh module
~~~~~~~~~~~~~~~~~~~~~~~~~

* **No in-place mutation of the mesh objects.**
  To make the maintenance of the mesh routines easier, all in-place transformation of the objects have been removed.
  Code such as the following::

    mesh = cpt.load_mesh("...")
    mesh.translate_x(1.0)
    mesh.keep_immersed_part()
    mesh.show()  # Show translated and clipped mesh

  can be rewritted as::

    mesh = cpt.load_mesh("...")
    mesh = mesh.translated_x(1.0)
    mesh = mesh.immersed_part()
    mesh.show()

  or more consisely::

    mesh = cpt.load_mesh("...").translated_x(1.0).immersed_part()
    mesh.show()

  In the latter version, each method call returns a new Python object of class ``Mesh`` which is the object now referred to by the variable name ``mesh``.
  The new version makes it less error prone to have complex workflows such as::

    full_mesh = cpt.load_mesh("...")
    full_mesh_draft = full_mesh.vertices[:, 2].min()
    for draft in [1.0, 2.0, 3.0]:
        mesh = full_mesh.translated_z(full_mesh_draft - draft).immersed_part()
    ...

  without any risk to overwriting the original ``full_mesh``.

  The only usage of in-place transformations is for performance critical part of the code.
  Given that most hydrodynamical meshes are usually below 100k faces, Capytaine's mesh class is usually not the performance bottleneck.
  Computation intensive mesh transformations should be done with a dedicated meshing tool and not directly in Capytaine anyway.
  If you find yourself nonetheless struggling with performance issues of the new mesh module in Capytaine, please open an issue on Github.

* **No in-place mutation of the body objects.**
  For the same reasons as above, the in-place transformations of the ``FloatingBody`` object have been removed.

  Code such as the following::

    body = cpt.FloatingBody(mesh=mesh)
    body.rotation_center = (0, 0, -1)
    body.add_all_rigid_body_dofs()
    body.keep_only_dofs(['Heave'])
    body.keep_immersed_part()
    body.translate([0, 1, 0])

  can be rewritten as::

    body = cpt.FloatingBody(
      mesh=mesh,
      dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -1), only=["Heave"])
    )
    body = body.immersed_part()
    body = body.translated([0, 1, 0])

* **Only a few built-in mesh loaders, but transparently use external libraries**
  Only the domain mesh file formats (Nemoh's, WAMIT's, Hydrostar's and HAMS's) are built-in in Capytaine.
  Loading a mesh in a general purpose file format such as GMSH or STL is still easy, assuming a third party library supporting this file format is installed (see :doc:`user_manual/mesh`).

* **No more mesh writers**
  To reduce the burden of maintenance, mesh writers have been removed.
  The mesh objects can be exported to external libraries that can write mesh files.
  As a consequence ``export_as_Nemoh_directory`` has been moved out of Capytaine.

* **Symmetries are only available around the main axis.**
  The ``Plane`` and ``Axis`` objects have been removed.
  Symmetric meshes using ``ReflectionSymmetricMesh`` can now only be defined for global symmetries across the ``'xOz'`` and ``'yOz'`` planes.
  Other planes and local symmetries used to be supported in the previous version of Capytaine, but they were making the implementation much more complicated for little practical gain, so it has been chosen for this new version to reduce the scope but make sure that this feature is well integrated with all the other features of Capytaine
  Similarly, transformations with the ``mirrored`` and ``rotated`` methods don't work with arbitrary plane or axis anymore.
  Most transformation can still be performed by combining translation, rotation and mirroring.
  More complex transformations should be done in a dedicated meshing software.

* Consequently to the removal of ``Plane`` and ``Axis``, the methods
  ``Mesh.clipped`` and ``FloatingBody.clipped`` don't take as argument a
  ``Plane`` object, but directly a point and a normal (as two 3-ple of floats).

* **Rotation symmetric meshes have been completely reworked.**
  They are now well integrated with other features such as irregular frequencies removal.
  The user interface has been changed since the experimental methods from previous versions.
  See :class:`~capytaine.meshes.symmetric_meshes.RotationSymmetricMesh`.
  The method to create a symmetric mesh from a profile of points has also changed, see :meth:`~capytaine.meshes.symmetric_meshes.RotationSymmetricMesh.from_profile_points`.
  Also the dihedral symmetry (nested reflection symmetry within a rotation symmetry) is now supported.

* **Prototype translation symmetry has been removed.**
  As a consequence, the ``translation_symmetry`` arguments of the mesh generations functions has been removed.

* Mesh quality checks have been reworked.

* If no name is provided, no generic name is given to the mesh, no name is used.
  Meshes' names are only useful to keep track of Python objects, since printing the full list of points and faces is not very convenient.

* Default 3D display now uses `pyvista <https://docs.pyvista.org/>`_ as a backend instead of raw ``vtk``. Please consider installing ``pyvista``.
  The ``matplotlib`` backend is also still available for static mesh viewing.
  Some keyword argument might have been changed to uniformize usage of the two 3D backends.
  Former animation tooling has been replaced with new tools using the `vedo <https://vedo.embl.es/>`_ library as a backend.
  They can be found in :mod:`~capytaine.ui.vedo_animations` with examples of usage in the cookbook.

* The barely-used and barely-documented ``geometric_center`` attribute of the mesh and the bodies have been removed.

* Support for optionally using quadratures from Quadpy has been removed.

Other major changes
~~~~~~~~~~~~~~~~~~~

* Add :class:`~capytaine.bodies.multibodies.Multibody` meant to represent a
  **multibody system**. For hydrodynamics, this is equivalent to the previous
  behavior of coalescing bodies together. For hydrostatics, the new class
  is slightly more powerful, for instance by being able to keep track of
  several center of buoyancy and center of mass. (:pull:`822`)

  Joining bodies with :meth:`~capytaine.bodies.bodies.FloatingBody.join_bodies`
  or ``+`` now creates a :class:`~capytaine.bodies.multibodies.Multibody` instance. It can be converted back to a
  :class:`~capytaine.bodies.bodies.FloatingBody` instance with::

    both = body_1 + body_2  # `both` is now a Multibody
    both = (body_1 + body_2).as_FloatingBody  # Recover former behavior of joining FloatingBody with a FloatingBody

  :class:`~capytaine.bodies.bodies.FloatingBody` and
  :class:`~capytaine.bodies.multibodies.Multibody` both inherits from the
  abstract class :class:`~capytaine.bodies.abstract_bodies.AbstractBody`
  (:pull:`873`)


* **Hydrostatics (including multibody hydrostatics) are now stored in an xarray Dataset**.
   The main way to get hydrostatics properties about a ``FloatingBody`` or a
   new ``Multibody`` is through :func:`~capytaine.compute_hydrostatics_dataset`.
   This function is called automatically by default when using
   :meth:`~capytaine.bem.solver.BEMSolver.fill_dataset`, such that hydrostatics
   data should always be found near hydrodynamics data in datasets outputted by
   Capytaine.
   The new data stored by default in the output dataset include notably the definition
   of the ``rotation_center`` of the bodies.
   For multibody cases, the dataset contains a new `body` dimension to store
   the rotation center, the center of mass and the center of buoyancy of each
   body individually (other magnitudes such as hydrostatic stiffness, inertia
   matrix, added mass, radiation damping and excitation forces are still stored
   using the ``radiating_dof`` and ``influenced_dof`` dimensions which may
   contains body names, but cannot be stored body-by-body in the ``body``
   dimension).
   (:pull:`887`)
   The xarray hydrostatics dataset can be expoted to Nemoh format with the new
   :func:`~capytaine.io.legacy.export_hydrostatics_from_dataset`, which is called automatically by :func:`~capytaine.io.xarray.export_dataset` when choosing Nemoh's format. (:pull:`898`)


* **New internal data model for rigid body dofs** with the classes
  :class:`~capytaine.bodies.dofs.TranslationDof` and
  :class:`~capytaine.bodies.dofs.RotationDof`. (:pull:`838`)
  This should ensure that rigid dofs are detected and treated as such, without
  relying only on their name.
  This is relevant for multiple bodies and articulated bodies when computing:
  * hydrostatics, where the exact hydrostatic stiffness formula for rigid body dofs can be used instead of the approximation for generalized dofs.
  * forward speed, where the m-term is currently only implemented for rigid dofs.

  Subsequently :func:`~capytaine.bodies.dofs.rigid_body_dofs` now instantiates the new dof class instead of a placeholder.
  It now has an additional input argument ``only`` to have only some of the six rigid body dofs.
  Also a rotation center should be passed, because the properties of the body
  (e.g. `center_of_mass`) cannot be accessed at that stage. (:pull:`838`)


* **Far-field mean drift force** evaluation with :func:`~capytaine.post_pro.mean_drift_force.far_field_mean_drift_force`.
  Export of far-field mean drift force to WAMIT `.8` format has also been implemented.

Minor changes
~~~~~~~~~~~~~

* Add option :code:`'lu_decomposition_with_overwrite'` for the :code:`linear_solver` of :class:`~capytaine.bem.engines.DefaultMatrixEngine`, which reduces the RAM usage of the solver (:pull:`775`).

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

* Mesh quadrature (if defined) is used for hydrostatics computations.
  Since meshes have no quadrature by default, the default hydrostatics results remain the same.
  Using `mesh.with_quadrature('Gauss-Legendre 2')` can lead to slightly better results,
  although the limitation is often the geometric approximation of the shape by a flat panels mesh.
  (:pull:`847`)

* **Breaking** The deprecated ``FreeSurface`` class and ``BEMSolver.get_potential_on_mesh`` and ``BEMSolver.get_free_surface_elevation`` methods have been removed.
  They can be replaced by just any mesh representation of the free surface and the methods :meth:`~capytaine.bem.solver.BEMSolver.compute_potential` and :meth:`~capytaine.bem.solvr.BEMSolver.compute_free_surface_elevation`.

* Computing the pressure or free surface elevation in post-processing does not
  allocate a very large matrix for a single matrix-vector product, but instead
  allocate and evaluate only a few rows of the matrix at a time (:pull:`860`).

* Warn the user if the ``water_depth`` would be better set to infinity than to a very large number (:pull:`880`).

* **Breaking** Remove the `body_name` dimension from output dataset, and the possibility to automatically stack together outputs for different "bodies" (that is geometries) in the same dataset.
  Now :func:`~capytaine.io.xarray.problems_from_dataset` and :meth:`~capytaine.bem.solver.fill_dataset` can only take a single body as input (or a list with a single element for backward compatibility).
  This decision is meant to make it easier to output multibody hydrostatics (see above), meshes data and local pressures in the output dataset. (:pull:`885`)

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
  computed displaced volume is always actually the volume below :math:`z=0` and the
  center of buoyancy is always below the free surface.
  The inertia matrix always uses the displaced water mass for the mass and
  compute the inertia moments on the full shape if the full mesh is provided.
  (:pull:`794`)

* Fix the frequency type in the dimensions of the dataset returned by :meth:`~cpt.io.xarray.kochin_data_array`. Previously, the dimension was always named ``omega``;
  it is now named ``omega``, ``freq``, ``period``,  ``wavenumber`` or ``wavelength`` depending on the user settings.

* Fix a bug when computing the Kochin function on a mesh with a lid. The Kochin function now also takes the faces on the lid into account. (:issue:`833`)

* Fix WAMIT output in ``.1`` file where 0 and infinite frequency added mass where interverted. (:pull:`855`)

* Fix resolution warning message to include ``lid_mesh`` (:issue:`867` and :pull:`868`)

* Progress bar display progress of groups of problems solved, not single problems, leading to more accurate time estimation (:pull:`882`)

* The RAO used to be computed wrongly with the transpose of the actual transfer function. (:issue:`891` and :pull:`892`)

* :func:`~capytaine.io.xarray.assemble_matrices` now returns the matrices ordered along the frequency axis like the input ``results``, instead of sorted by increasing frequency. (:issue:`797`)

* Exporting to Nemoh's ``Hydrostatics.dat`` stores the displaced volume, not the total volume of the body.

* The logging is not forcefully overridden by Capytaine on import. Capytaine
  logging's setting are only set if no logging has been set before.
  (:pull:`922`)

Internals
~~~~~~~~~

* **Breaking** The ``green_function`` is not an attribute of the :class:`~capytaine.bem.solver.BEMSolver` anymore, but of the engine.
  The motivation is that not all engines can be made compatible with all Green function implementations (although the builtins one are).
  The possibility to call ``BEMSolver(green_function=...)`` is kept as a convenient shortcut to ``BEMSolver(engine=DefaultMatrixEngine(green_function=...))``.
  Calls to ``BEMSolver(green_function=..., engine=...)`` now raise an error. (:pull:`752`)
  Post-processing new requires the implementation of the methods ``build_S_matrix`` and ``build_fullK_matrix`` by the engine (:pull:`753`)

* New implementation of the block symmetric matrices for mesh symmetry, now
  used by :class:`~capytaine.bem.engines.DefaultMatrixEngine` (:pull:`754`).

* Rafactor of the :class:`~capytaine.bem.engines.DefaultMatrixEngine` to make the
  caching more straightforward and improve its interaction with LU
  decomposition and symmetries. (:pull:`755`)

* The whole ``matrices`` module as well as the corresponding engines in
  ``bem/engines.py`` and ``tool/lru_caches.py`` have been removed from
  Capytaine. For compatibility, they will remain accessible from a separate
  package. (:pull:`757`, :pull:`765`)

* The new ``Mesh`` class has a ``faces_metadata`` attribute storing fields
  defined on each faces of the mesh. When faces are added or removed, the
  metadata are automatically updated accordingly. (:pull:`791`)

* Move hydrostatics routines in a dedicated module and rewrite corresponding tests (:pull:`794`)
  Change behavior of hydrostatic stiffness when no dof is defined: return an empty matrix rather than raising an error, in order to better fit in multibody cases (:pull:`883`).

* Refactor the implementation of the timer to make it easier to include more steps (:pull:`809`)

* Parameters ``free_surface``, ``water_depth`` and ``wavenumber`` are always
  keyword arguments in the engine and the Green function (:pull:`812`)

* **Breaking** :meth:`~capytaine.bodies.FloatingBody.add_rotation_dof` takes
  arguments called ``rotation_center`` and ``direction`` instead of an ``Axis``
  object. Also the geometric center is not used anymore as a fallback value for
  ``rotation_center``.

* The S matrix can now be computed even if the K matrix is not defined,
  for example when evaluating the :meth:`~capytaine.bem.solver.compute_free_surface_elevation` along the waterline.

* New function implemented :func:`~capytaine.meshes.waterline_integral` to compute the integral of a function along the water line of a mesh.

* **Breaking** "Full double-layer matrices" are now stored as array of shape (3, n, m)
  (or 3 matrices of shape (n, m)) instead of (n, m, 3) (:pull:`869`).
