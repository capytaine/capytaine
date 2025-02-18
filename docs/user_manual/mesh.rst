======
Meshes
======

Naming
------

Meshes in Capytaine have a ``name`` attribute.
It is optional and is mostly used for clearer logging and outputs.
A ``name`` optional argument can be provided to all methods below to initialize
a mesh or transform a mesh to set the name of the new mesh.


.. note::
   A mesh in Capytaine is merely a set of independent faces (triangles or quadrangles).
   Connectivities are not required for the resolution.
   Having a mesh that is not watertight, with small gaps between the faces or a
   few missing faces, does not lead to qualitatively different results.


Initialization
--------------

Importing with included Meshmagick readers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To load an existing mesh file, use the following syntax::

    import capytaine as cpt

    mesh = cpt.load_mesh('path/to/mesh.dat', file_format='nemoh')
    body = cpt.FloatingBody(mesh=mesh)

The above example uses Nemoh's mesh format, as defined e.g. on page 19 of `Nemoh v3.0.1 manual`_.

.. _`Nemoh v3.0.1 manual`: https://gitlab.com/api/v4/projects/41313230/packages/generic/nemoh/v3.0.1/Nemoh_Manual_v3.0.1.pdf

Thanks to inclusion of code from `Meshmagick <https://github.com/lheea/meshmagick/>`_,
numerous other mesh format can be imported.
The file format can be given with the :code:`file_format` optional argument.
If no format is given, the code will try to infer it from the file extension::

    mesh = cpt.load_mesh('path/to/mesh.msh')  # gmsh file

The formats currently supported in reading are listed in the following table (adapted from the documentation of Meshmagick).

+-----------+-----------------+----------------------+-----------------+
| File      | Software        | Keywords             | Extra features  |
| extension |                 |                      |                 |
+===========+=================+======================+=================+
|   .mar    | NEMOH [#f1]_    | nemoh, mar           | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .nem    | NEMOH [#f1]_    | nemoh_mesh, nem      |                 |
+-----------+-----------------+----------------------+-----------------+
|   .gdf    | WAMIT [#f2]_    | wamit, gdf           | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .inp    | DIODORE [#f3]_  | diodore-inp, inp     |                 |
+-----------+-----------------+----------------------+-----------------+
|   .DAT    | DIODORE [#f3]_  | diodore-dat          |                 |
+-----------+-----------------+----------------------+-----------------+
|   .pnl    | HAMS            | pnl, hams            | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .hst    | HYDROSTAR [#f4]_| hydrostar, hst       | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .nat    |    -            | natural, nat         |                 |
+-----------+-----------------+----------------------+-----------------+
|   .msh    | GMSH [#f5]_     | gmsh, msh            |                 |
+-----------+-----------------+----------------------+-----------------+
|   .rad    | RADIOSS         | rad, radioss         |                 |
+-----------+-----------------+----------------------+-----------------+
|   .stl    |    -            | stl                  |                 |
+-----------+-----------------+----------------------+-----------------+
|   .vtu    | PARAVIEW [#f6]_ | vtu                  |                 |
+-----------+-----------------+----------------------+-----------------+
|   .vtp    | PARAVIEW [#f6]_ | vtp                  |                 |
+-----------+-----------------+----------------------+-----------------+
|   .vtk    | PARAVIEW [#f6]_ | paraview-legacy, vtk |                 |
+-----------+-----------------+----------------------+-----------------+
|   .tec    | TECPLOT [#f7]_  | tecplot, tec         |                 |
+-----------+-----------------+----------------------+-----------------+
|   .med    | SALOME [#f8]_   | med, salome          |                 |
+-----------+-----------------+----------------------+-----------------+

.. [#f1] NEMOH is an open source BEM Software for seakeeping developed at
         Ecole Centrale de Nantes (LHEEA)
.. [#f2] WAMIT is a BEM Software for seakeeping developed by WAMIT, Inc.
.. [#f3] DIODORE is a BEM Software for seakeeping developed by PRINCIPIA
.. [#f4] HYDROSTAR is a BEM Software for seakeeping developed by
         BUREAU VERITAS
.. [#f5] GMSH is an open source meshing software developed by C. Geuzaine
         and J.-F. Remacle. Version 4 of the file format requires meshio
         be installed independently.
.. [#f6] PARAVIEW is an open source visualization software developed by
         Kitware
.. [#f7] TECPLOT is a visualization software developed by Tecplot
.. [#f8] SALOME-MECA is an open source software for computational mechanics
         developed by EDF-R&D


Not all metadata is taken into account when reading the mesh file.
For instance, the body symmetry is taken into account only for the ``.mar``, ``.pnl``, ``.gdf`` and ``.hst`` file formats.
Feel free to open an issue on Github to suggest improvements.


Importing with Meshio
~~~~~~~~~~~~~~~~~~~~~

Mesh can also be imported using the `meshio <https://pypi.org/project/meshio/>`_
library. Unlike the Meshmagick mesh readers mentioned above, this library is
not packaged with Capytaine and need to be installed independently::

    pip install meshio

The `meshio` mesh object can converted to Capytaine's mesh
format with the :func:`~capytaine.io.meshio.load_from_meshio` function::

    from capytaine.io.meshio import load_from_meshio
    cpt_mesh = cpt.load_from_meshio(mesh, name="My mesh")

This features allows to use `pygmsh <https://pypi.org/project/pygmsh/>`_ to
generate the mesh, since this library returns mesh in the same format as meshio.
Below is an example of a mesh generation with `pygmsh` (which also needs to be
installed independently)::

    import pygmsh
    offset = 1e-2
    T1 = 0.16
    T2 = 0.37
    r1 = 0.88
    r2 = 0.35
    with pygmsh.occ.Geometry() as geom:
        cyl = geom.add_cylinder([0, 0, 0], [0, 0, -T1],  r1)
        cone = geom.add_cone([0, 0, -T1], [0, 0, -T2], r1, r2)
        geom.translate(cyl, [0, 0, offset])
        geom.translate(cone, [0, 0, offset])
        geom.boolean_union([cyl, cone])
        gmsh_mesh = geom.generate_mesh(dim=2)
    mesh = cpt.load_mesh(gmsh_mesh, name="my_pygmsh_mesh")


Predefined simple shapes
~~~~~~~~~~~~~~~~~~~~~~~~

Capytaine include mesh generators for a few simple shapes. They are mostly
meant for teaching (they are extensively used in the examples of this
documentation) as well as for testing.
The most useful ones are
:func:`~capytaine.meshes.predefined.spheres.mesh_sphere`,
:func:`~capytaine.meshes.predefined.cylinders.mesh_vertical_cylinder`,
:func:`~capytaine.meshes.predefined.cylinders.mesh_horizontal_cylinder`,
:func:`~capytaine.meshes.predefined.rectangles.mesh_parallelepiped`.
Some applications may also make use of flat shapes
:func:`~capytaine.meshes.predefined.cylinders.mesh_disk` and
:func:`~capytaine.meshes.predefined.rectangles.mesh_rectangle`.
Refer to their documentation for details about the parameters they accepts.

Since version 2.1, their resolution can be set by the ``faces_max_radius``
parameter which specifies the maximal size of a face in the mesh.

.. note::
    There are several ways to measure the size of a face and the resolution of a mesh.
    In Capytaine, the size of faces is usually quantified with the *radius* of the face, that is the maximal distance between the center of the face and its vertices.
    The resolution of a mesh is estimated as the maximal radius among all the faces in the mesh, that is the radius of the biggest face.


Creating from scratch
~~~~~~~~~~~~~~~~~~~~~

Alternatively, a mesh can be defined by giving a list of vertices and faces::

    mesh = cpt.Mesh(vertices=..., faces=..., name="my_mesh")

The vertices are expected to be provided as a Numpy array of floats with shape ``(nb_vertices, 3)``.
The faces are provided as a Numpy array of ints with shape ``(nb_faces, 4)``, such that the four integers on a line are the indices of the vertices composing that face::

    v = np.array([[0.0, 0.0, -1.0],
                  [1.0, 0.0, -1.0],
                  [1.0, 1.0, -1.0],
                  [0.0, 1.0, -1.0]])
    f = np.array([[0, 1, 2, 3]])
    single_face_mesh = cpt.Mesh(vertices=v, faces=f)

The ordering of the vertices define the direction of the normal vector, using
normal right rotation. In other words, the normal vector is towards you if you
see the vertices as being in counterclockwise order. In the above example, the
normal vector is going up.

Triangular faces are supported as quadrilateral faces with the same vertex
repeated twice::

    single_triangle_mesh = cpt.Mesh(vertices=v, faces=np.array([[0, 1, 2, 2]]))


Creating a symmetric mesh
~~~~~~~~~~~~~~~~~~~~~~~~~

Several mesh symmetries can be used by Capytaine to speed up the computation.
The most useful one is the vertical plane symmetry.
A mesh with such a symmetry is stored by Capytaine with the
:class:`~capytaine.meshes.symmetric.ReflectionSymmetricMesh` class.
It is defined with an other mesh of the half and a plane (and optionally a name
like the usual meshes)::

    half_mesh = cpt.load_mesh(...)
    mesh = cpt.ReflectionSymmetricMesh(half_mesh, cpt.xOz_Plane, name="my full mesh")

Two vertical plane symmetries can be nested to be used by Capytaine (assuming
that the two planes are orthogonal)::

    quarter_mesh = cpt.load_mesh(...)
    half_mesh = cpt.ReflectionSymmetricMesh(half_mesh, cpt.yOz_Plane)
    mesh = cpt.ReflectionSymmetricMesh(half_mesh, cpt.xOz_Plane)

All the method defined afterwards in this documentation should be applicable
for ``ReflectionSymmetricMesh`` as well as for standard ``Mesh``.

You can consider using the ``clipped`` method discussed below to create a symmetric mesh::

    half_mesh = original_mesh.clipped(plane=cpt.xOz_Plane)
    mesh = cpt.ReflectionSymmetricMesh(half_mesh, cpt.xOz_Plane)


Display
-------

Use the :code:`show` method to display the mesh in 3D using VTK (if installed)
with the :meth:`~capytaine.meshes.meshes.Mesh.show`::

    mesh.show()

or with Matplotlib (if installed) with
:meth:`~capytaine.meshes.meshes.Mesh.show_matplotlib`::

    mesh.show_matplotlib()


Geometric transformations
-------------------------
Several functions are available to transform existing meshes.

Below is a list of most of the available methods.
All of them can be applied to both meshes or to floating bodies, in which case
the degrees of freedom will also be transformed::

    # TRANSLATIONS
    mesh.translated_x(10.0)
    mesh.translated_y(10.0)
    mesh.translated_z(10.0)
    mesh.translated([10.0, 5.0, 2.0])

    # Translation such that point_a would become equal to point_b
    mesh.translated_point_to_point(point_a=[5, 6, 7], point_b=[4, 3, 2])

    # ROTATIONS
    mesh.rotated_x(3.14/5)  # Rotation of pi/5 around the Ox axis
    mesh.rotated_y(3.14/5)  # Rotation of pi/5 around the Oy axis
    mesh.rotated_z(3.14/5)  # Rotation of pi/5 around the Oz axis

    # Rotation of pi/5 around an arbitrary axis.
    from capytaine import Axis
    my_axis = Axis(vector=[1, 1, 1], point=[3, 4, 5])
    mesh.rotated(axis=my_axis, angle=3.14/5)

    # Rotation around a point such that vec1 would become equal to vec2
    mesh.rotated_around_center_to_align_vector(
        center=(0, 0, 0),
        vec1=(1, 4, 7),
        vec2=(9, 2, 1)
    )

    # REFLECTIONS
    from capytaine import Plane
    mesh.mirrored(Plane(normal=[1, 2, 1], point=[0, 4, 5]))

All the above methods can also be applied to :class:`~capytaine.meshes.geometry.Plane`
and :class:`~capytaine.meshes.geometry.Axis` objects.

Meshes can also be merged together with the :code:`+` operator::

    larger_mesh = mesh_1 + mesh_2

Finally, meshes can be clipped with a :class:`~capytaine.meshes.geometry.Plane`.
The plane is defined by a point belonging to it and a normal vector::

    xOy_Plane = Plane(point=(0, 0, 0), normal=(0, 0, 1))
    clipped_mesh = mesh.clipped(xOy_Plane)

Beware that the orientation of the normal vector of the :code:`Plane` will
determine which part of the mesh will be returned::

    higher_part = mesh.clipped(Plane(point=(0, 0, 0), normal=(0, 0, -1)))
    lower_part = mesh.clipped(Plane(point=(0, 0, 0), normal=(0, 0, 1)))
    # mesh = lower_part + higher_part

The method :meth:`immersed_part` will clip the body with respect to two
horizontal planes at :math:`z=0` and :math:`z=-h`::

    clipped_body = mesh.immersed_part(water_depth=10)

.. note::
    Most transformation methods exist in two versions:

    * one, named as a infinitive verb (`translate`, `rotate`, `clip`,
      `keep_immersed_part`, ...), is an in-place transformation;
    * the other, named as a past participle (`translated`, `rotated`,
      `clipped`, `immersed_part`, ...), is the same transformation but
      returning a new object.

    In most cases, performance is not significant and the method returning a
    new object should be preferred. In-place transformation are currently kept
    for backward compatibility, but they make the code significantly more
    complicated and their removal might be considered in the future.


Extracting or generating a lid
------------------------------

If you loaded a mesh file already containing a lid on the :math:`z=0` plane,
the hull and the lid can be split with the
:meth:`~capytaine.meshes.meshes.Mesh.extract_lid` method::

    full_mesh = cpt.load_mesh(...)
    hull_mesh, lid_mesh = full_mesh.extract_lid()

If your mesh does not have a lid, and you'd like to have irregular frequencies
removal, you can generate a lid using
:meth:`~capytaine.meshes.meshes.Mesh.generate_lid` as follows::

    lid_mesh = hull_mesh.generate_lid()

The mesh is generated on the free surface by default.
Since support for panels on the free surface is still experimental, it might be
more robust (but less efficient) to define a lid slightly below the free surface::

    lid_mesh = hull_mesh.generate_lid(z=-0.1)

The lower the lid, the more robust the computation, but also the less
irregular frequencies are removed. The method
:meth:`~capytaine.meshes.meshes.Mesh.lowest_lid_position` estimates the lowest
position such that all irregular frequencies below a given frequency are removed::

    lid_mesh = hull_mesh.generate_lid(z=hull_mesh.lowest_lid_position(omega_max=10.0))

The method :meth:`~capytaine.meshes.meshes.Mesh.extract_lid` also accepts an
optional argument ``faces_max_radius`` to set the resolution of the lid. By
default, the mean resolution of the hull mesh is used.

See :doc:`body` for detail on how to assign a lid mesh when defining a floating
body.

.. note::
   The lid does not need neither to cover the whole interior free surface, nor
   to be connected with the hull mesh. The lid automatically generated by
   :meth:`~capytaine.meshes.meshes.Mesh.extract_lid` typically does not.
   Nonetheless, the more interior free surface is covered, the more efficiently
   the irregular frequencies will be removed.



Defining an integration quadrature
----------------------------------

.. note::
   Quadratures are an advanced feature meant to experiment with numerical schemes.
   The best compromise between precision and performance is often not to bother
   with it and keep the default integration scheme.

During the resolution of the BEM problem, the Green function has to be
integrated on each panel of the mesh. Parts of the Green function (such as the
:math:`1/r` Rankine terms) are integrated using an exact analytical expression
for the integral. Other parts of the Green function rely on numerical
integration. By default, this numerical integration is done by taking the value
at the center of the panel and multiplying by its area. For a more accurate
intagration, an higher order quadrature can be defined.

To define a quadrature scheme for a mesh, run the following command::

    mesh.compute_quadrature(method="Gauss-Legendre 2")

The quadrature data can then be accessed at::

    mesh.quadrature_points

and will be used automatically when needed.

.. warning:: Transformations of the mesh (merging, clipping, ...) may reset the quadrature.
             Compute it only on your final mesh.

.. warning:: Quadratures schemes have been designed with quadrilateral panels.
             They work on triangular panels, but might not be as optimal then.

Alternatively, the :func:`~capytaine.meshes.meshes.Mesh.compute_quadrature`
also accepts methods from the `Quadpy` package::

    import quadpy
    mesh.compute_quadrature(method=quadpy.c2.get_good_scheme(8))
