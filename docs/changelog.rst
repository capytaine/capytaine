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
New in version 3.0 (2026-??-??)
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

  In the latter version, each line returns a new Python object of class ``Mesh`` which is the object now referred to by the variable name ``mesh``.
  The new version makes it less error prone to have complex workflows such as::

    full_mesh = cpt.load_mesh("...")
    full_mesh_draft = full_mesh.vertices[:, 2].min()
    for draft in [1.0, 2.0, 3.0]:
        mesh = full_mesh.translated_z(full_mesh_draft - draft).immersed_part()
    ...

  without any risk to overwriting the original ``full_mesh``.

  Subsequently, the ``copy()`` method has been removed.

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
      dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -1))
    )
    body = body.with_only_dofs(['Heave'])
    body = body.immersed_part()
    body = body.translated([0, 1, 0])

* **Only a few built-in mesh loaders, but transparently use external libraries**
  Only the domain mesh file formats (Nemoh's, WAMIT's, Hydrostar's and HAMS's) are built-in in Capytaine.
  Loading a mesh in a general purpose file format such as GMSH or STL is still easy, assuming a third party library supporting this file format is installed (see :doc:`user_manual/mesh`).

* **No more mesh writers**
  To reduce the burden of maintenance, mesh writers have been removed, but the
  mesh objects can be exported to external libraries that can write mesh files.
  As a consequence ``export_as_Nemoh_directory`` has been moved out of Capytaine.

* **Symmetries are only available around the main axis.**
  The ``Plane`` and ``Axis`` objects have been removed.
  Symmetric meshes using ``ReflectionSymmetricMesh`` can now only be defined for global symmetries across the ``'xOz'`` and ``'yOz'`` planes.
  Other planes and local symmetries used to be supported in the previous version of Capytaine, but they were making the implementation much more complicated for little practical gain, so it has been chosen for this new version to reduce the scope but make sure that this feature is well integrated with all the other features of Capytaine
  Similarly, transformations with the ``mirrored`` and ``rotated`` methods don't work with arbitrary plane or axis anymore.
  Most transformation can still be performed by combining translation, rotation and mirroring.
  More complex transformations should be done in a dedicated meshing software.

* **Rotation symmetric meshes have been completely reworked.**
  They are now well integrated with the other features such as irregular frequencies removal.
  The user interface has been changed since the experimental methods from previous versions.
  See :class:`~capytaine.meshes.symmetric_meshes.RotationSymmetricMesh`.
  The method to create a symmetric mesh from a profile of points has also changed, see :meth:`~capytaine.meshes.symmetric_meshes.RotationSymmetricMesh.from_profile_points`.
  Also the dihedral symmetry (nested reflection symmetry within a rotation symmetry) is now supported.

* **Prototype translation symmetry has been removed.**
  As a consequence, the ``translation_symmetry`` arguments of the mesh generations functions has been removed.

* Different quality checks suite for given meshes.

* ``Mesh.clipped`` and ``FloatingBody.clipped`` don't take as argument a
  ``Plane`` object, but directly a point and a normal (as two 3-ple of floats).

* If no name is provided, no generic name is given to the mesh, no name is used.
  Meshes' names are only useful to keep track of Python objects, since printing the full list of points and faces is not very convenient.

* Default 3D display now uses `pyvista <https://docs.pyvista.org/>`_ as a backend instead of raw ``vtk``. Please consider installing ``pyvista``.
  The ``matplotlib`` backend is also still available for static mesh viewing.
  Some keyword argument might have been changed to uniformize usage of the two 3D backends.
  Former animation tooling has been replaced with new tools using the `vedo <https://vedo.embl.es/>`_ library as a backend. They can be found in :mod:`~capytaine.ui.vedo_animations` with examples of usage in the cookbook.

* The barely-used and barely-documented ``geometric_center`` attribute of the mesh and the bodies have been removed.

* Support for optionally using quadratures from Quadpy has been removed.

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

  :class:`~capytaine.bodies.bodies.FloatingBody` and
  :class:`~capytaine.bodies.multibodies.Multibody` both inherits from the
  abstract class :class:`~capytaine.bodies.abstract_bodies.AbstractBody`
  (:pull:`873`)

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

* Mesh quadrature (if defined) is used for hydrostatics computations.
  Since meshes have no quadrature by default, the default hydrostatics results remain the same.
  Using `mesh.with_quadrature('Gauss-Legendre 2')` can lead to slightly better results,
  although the limitation is often the geometric approximation of the shape by a flat panels mesh.
  (:pull:`847`)
  
* Add function :func:`~capytaine.post_pro.mean_drift_force.far_field_mean_drift_force` to compute the horizontal mean drift forces using far field formulation.
  Only the single-direction second-order term is currently implemented.

* **Breaking** The deprecated ``FreeSurface`` class and ``BEMSolver.get_potential_on_mesh`` and ``BEMSolver.get_free_surface_elevation`` methods have been removed.
  They can be replaced by just any mesh representation of the free surface and the methods :meth:`~capytaine.bem.solver.BEMSolver.compute_potential` and :meth:`~capytaine.bem.solvr.BEMSolver.compute_free_surface_elevation`.

* Computing the pressure or free surface elevation in post-processing does not allocate a very large matrix for a single matrix-vector product, but instead allocate and evaluate only a few rows of the matrix at a time (:pull:`860`).


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

* New function implemented :func:`~capytaine.meshes.waterline_integral` to compute the integral of a function along the water line of a mesh.

* "Full double-layer matrices" are now stored as array of shape (3, n, m)
  (or 3 matrices of shape (n, m)) instead of (n, m, 3) (:pull:`869`).
