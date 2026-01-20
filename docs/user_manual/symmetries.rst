Mesh symmetries
===============

Defining a symmetric mesh
~~~~~~~~~~~~~~~~~~~~~~~~~

Mesh symmetries can be used to speed up the computation.
Four kind of symmetries are supported by Capytaine:

* **Single plane symmetry** with respect to the :math:`x = 0` or :math:`y = 0` plane.
  This is the symmetry of most ship hulls and is thus implemented in almost all linear sea-keeping codes.

    A mesh with such a symmetry is stored by Capytaine with the
    :class:`~capytaine.new_meshes.symmetric_meshes.ReflectionSymmetricMesh` class.
    It is defined with an other mesh of the half and a plane (and optionally a name
    like the usual meshes)::

        half_mesh = cpt.load_mesh(...)
        mesh = cpt.ReflectionSymmetricMesh(half_mesh, plane="xOz", name="my full mesh")

    When loading a file in a format that support defining a symmetry (`gdf`,
    `hst`, `mar`, `pnl`), the ``ReflectionSymmetricMesh`` is returned
    automatically by ``load_mesh``.

* **Two plane symmetries** with respect to both the :math:`x = 0` and :math:`y = 0` plane.

    Two vertical plane symmetries can be nested to be used by Capytaine (assuming
    that the two planes are orthogonal)::

        quarter_mesh = cpt.load_mesh(...)
        half_mesh = cpt.ReflectionSymmetricMesh(half_mesh, plane="yOz")
        mesh = cpt.ReflectionSymmetricMesh(half_mesh, plane="x0z")

* **Rotation symmetry** of a shape repeated by rotation around the :math:`z`-axis any number of times.

    It can be defined either from the repetition of an existing mesh::

        wedge_mesh = cpt.load_mesh(...)
        mesh = cpt.RotationSymmetricMesh(wedge=wedge_mesh, n=4)

    or for an axysymmetric geometry from a list of points along a "meridian" of the shape::

        meridian_points = np.array([(np.sqrt(1-z**2), 0.0, z) for z in np.linspace(-1.0, 1.0, 10)])
        sphere = cpt.RotationSymmetricMesh.from_profile_points(meridian_points, n=10)

..
    * **Dihedral symmetry** is the combination of a plane symmetry and the rotation symmetry.

        It is defined by nesting a ``RotationSymmetricMesh`` into a ``ReflectionSymmetricMesh``::

            half_wedge = cpt.load_mesh(...)
            inner_mesh = cpt.RotationSymmetricMesh(half_wedge, n=4)
            mesh = cpt.ReflectionSymmetricMesh(half=inner_mesh, plane="xOz")

        The nesting in the other order is supported, but not as efficient.

Manipulating a symmetric mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All the methods defined in the :doc:`mesh` section can be used on the ``ReflectionSymmetricMesh`` and ``RotationSymmetricMesh``.

It the resulting object is not symmetric anymore, the symmetry is lost and a bar ``Mesh`` of the whole surface is returned::

    mesh = cpt.ReflectionSymmetricMesh(..., plane="xOz")
    x_translated_mesh = mesh.translated_x(1.0)
    print(x_translated_mesh.__class__)  # ReflectionSymmetricMesh
    y_translated_mesh = mesh.translated_y(1.0)
    print(y_translated_mesh.__class__)  # Mesh

In particular, joining meshes with ``+`` or ``join_meshes`` conserves the symmetry assuming all meshes have the same symmetry.
Clipping the part of the mesh above the free surface with ``immersed_part`` should always conserve the symmetry.


Using a symmetric mesh to define a floating body
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The symmetric mesh can be used to setup a floating body::

    mesh = cpt.ReflectionSymmetricMesh(..., plane="xOz")
    lid_mesh = cpt.ReflectionSymmetricMesh(..., plane="xOz")
    body = cpt.FloatingBody(
                mesh=mesh,
                lid_mesh=lid_mesh,
                dofs=cpt.rigid_body_dofs()
                )

.. warning::
   When using a lid mesh for irregular frequencies removal, the lid mesh and
   the hull mesh should have the same symmetry, otherwise the symmetry will be
   ignored when solving the BEM problem.
   The methods ``generate_lid`` and ``extract_lid`` preserve the symmetry,
   although the generated lid might not be optimal when the stub of the
   symmetric mesh is small.

.. note::
   For all the symmetries described above, the **mesh** (and the lid mesh)
   needs to be symmetric, but the full problem might not be.
   In other word, even when working with a symmetric mesh, it is possible to
   consider incoming waves from any directions, or generalized degrees of
   freedom that are not symmetric.


Expected performance gain
~~~~~~~~~~~~~~~~~~~~~~~~~

Asymptotically for large problems:

+---------------------+-----------------+-----------+
| Symmetry            | Resolution time | RAM usage |
+---------------------+-----------------+-----------+
| 1 plane symmetry    | 1/2             | 1/2       |
| 2 plane symmetries  | 1/4             | 1/4       |
| n rotation symmetry | 1/n             | 1/n       |
+---------------------+-----------------+-----------+

Note that the theoretical performance gain described above might not be reached
in practice, especially for smaller problems.
For instance, the threading parallelisation is currently less efficient on
highly symmetric meshes.
