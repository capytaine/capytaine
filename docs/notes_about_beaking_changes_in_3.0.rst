
New meshes module
-----------------

Main differences
~~~~~~~~~~~~~~~~

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

* **Symmetries are only available around the main axis.**
  The ``Plane`` and ``Axis`` objects have been removed.
  Symmetric meshes using ``ReflectionSymmetricMesh`` can now only be defined for global symmetries across the ``'xOz'`` and ``'yOz'`` planes.
  Other planes and local symmetries used to be supported in the previous version of Capytaine, but they were making the implementation much more complicated for little practical gain, so it has been chosen for this new version to reduce the scope but make sure that this feature is well integrated with all the other features of Capytaine
  Similarly, transformations with the ``mirrored`` and ``rotated`` methods don't work with arbitrary plane or axis anymore.
  Most transformation can still be performed by combining translation, rotation and mirroring.
  More complex transformations should be done in a dedicated meshing software.

* Different quality checks suite for given meshes.

* ``Mesh.clipped`` and ``FloatingBody.clipped`` don't take as argument a
  ``Plane`` object, but directly a point and a normal (as two 3-ple of floats).

* If no name is provided, no generic name is given to the mesh, no name is used.
  Meshes' names are only useful to keep track of Python objects, since printing the full list of points and faces is not very convenient.

* Default 3D display now uses ``pyvista`` as a backend instead of raw ``vtk``. Please consider installing ``pyvista``.
  Animation support has not been implemented in this new backend yet, only static mesh viewing is available.
  The ``matplotlib`` backend is also still available for static mesh viewing.
  Some keyword argument might have been changed to uniformize usage of the two 3D backends.
