==========================
Meshes and floating bodies
==========================

.. contents:: Content

Importing a mesh with Meshmagick
--------------------------------

To create a new body using an existing mesh file, use the following syntax::

    from capytaine import FloatingBody

    body = FloatingBody.from_file('path/to/mesh.dat', file_format='nemoh')

The above example uses `Nemoh's mesh format`_.

.. _`Nemoh's mesh format`: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp

Thanks to Meshmagick, numerous other mesh format can be imported.
The file format can be given with the :code:`file_format` optional argument.
If no format is given, the code will try to infer it from the file extension::

    body = FloatingBody.from_file('path/to/mesh.msh')  # gmsh file

The formats currently supported in reading are listed in the following table (adapted from the documentation of Meshmagick).

+-----------+-----------------+----------------------+-----------------+
| File      | Software        | Keywords             | Extra features  |
| extension |                 |                      |                 |
+===========+=================+======================+=================+
|   .mar    | NEMOH [#f1]_    | nemoh, mar           | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .nem    | NEMOH [#f1]_    | nemoh_mesh, nem      |                 |
+-----------+-----------------+----------------------+-----------------+
|   .gdf    | WAMIT [#f2]_    | wamit, gdf           |                 |
+-----------+-----------------+----------------------+-----------------+
|   .inp    | DIODORE [#f3]_  | diodore-inp, inp     |                 |
+-----------+-----------------+----------------------+-----------------+
|   .DAT    | DIODORE [#f3]_  | diodore-dat          |                 |
+-----------+-----------------+----------------------+-----------------+
|   .hst    | HYDROSTAR [#f4]_| hydrostar, hst       | Symmetries      |
+-----------+-----------------+----------------------+-----------------+
|   .nat    |    -            | natural, nat         |                 |
+-----------+-----------------+----------------------+-----------------+
|   .msh    | GMSH 2 [#f5]_   | gmsh, msh            |                 |
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
         and J.-F. Remacle. Version 4 of the file format is not supported at the
         moment.
.. [#f6] PARAVIEW is an open source visualization software developed by
         Kitware
.. [#f7] TECPLOT is a visualization software developed by Tecplot
.. [#f8] SALOME-MECA is an open source software for computational mechanics
         developed by EDF-R&D


Not all metadata is taken into account when reading the mesh file.
For instance, the body symmetry is taken into account only for the `.mar` and `.hst` file formats.
Feel free to open an issue on Github to suggest improvements.


Importing a mesh with Meshio
----------------------------

Mesh can also be imported using the `meshio <https://pypi.org/project/meshio/>`_
library. Unlike the Meshmagick mesh readers mentioned above, this library is
not packaged with Capytaine and need to be installed independently::

    pip install meshio

A `meshio` mesh object can be read using the :code:`FloatingBody.from_meshio`
method::

    import meshio
    mesh = meshio.read("myfile.stl")
    body = FloatingBody.from_meshio(mesh, name="My floating body")

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
        mesh = geom.generate_mesh(dim=2)
    body = FloatingBody.from_meshio(mesh)


Display and animation
---------------------
Use the :code:`show` method to display the body in 3D using VTK (if installed)::

    body.show()

The :code:`animate` method can be used to visualize a given motion of the body::

    anim = body.animate(motion={"Heave": 0.1, "Surge": 0.1j}, loop_duration=1.0)
    anim.run()

The above example will present an interactive animation of the linear combination of heave and surge.

Jupyter notebooks can also include a (non-interactive) video of the animation::

    anim.embed_in_notebook(camera_position=(-1.0, -1.0, 1.0), resolution=(400, 300))


Geometric transformations
-------------------------
Several functions are available to transform existing bodies and meshes.

Most transformation methods exist in two versions: 

* one, named as a infinitive verb (`translate`, `rotate`, ...), is an in-place transformation;
* the other, named as a past participle (`translated`, `rotated`, ...), is the
  same transformation but returning a new object. 

In most cases, performance is not significant and the latter method should be
preferred since it makes code slightly easier to debug.

Below is a list of most of the available methods.
All of them can be applied to both meshes or to floating bodies, in which case
the degrees of freedom will also be transformed::

    # TRANSLATIONS
    body.translated_x(10.0)
    body.translated_y(10.0)
    body.translated_z(10.0)
    body.translated([10.0, 5.0, 2.0])

    # Translation such that point_a would become equal to point_b
    body.translated_point_to_point(point_a=[5, 6, 7], point_b=[4, 3, 2])

    # ROTATIONS
    body.rotated_x(3.14/5)  # Rotation of pi/5 around the Ox axis
    body.rotated_y(3.14/5)  # Rotation of pi/5 around the Oy axis
    body.rotated_z(3.14/5)  # Rotation of pi/5 around the Oz axis

    # Rotation of pi/5 around an arbitrary axis.
    from capytaine import Axis
    my_axis = Axis(vector=[1, 1, 1], point=[3, 4, 5])
    body.rotated(axis=my_axis, angle=3.14/5)

    # Rotation around a point such that vec1 would become equal to vec2
    body.rotated_around_center_to_align_vector(
        center=(0, 0, 0),
        vec1=(1, 4, 7),
        vec2=(9, 2, 1)
    )

    # REFLECTIONS
    from capytaine import Plane
    body.mirrored(Plane(normal=[1, 2, 1], point=[0, 4, 5]))

All the above method can also be applied to :class:`~capytaine.meshes.geometry.Plane`
and :class:`~capytaine.meshes.geometry.Axis` objects.


Joining
-------
Meshes and bodies can be merged together with the :code:`+` operator::

    both_bodies = body_1 + body_2

The :code:`+` operation is associative, that is :code:`(body_1 + body_2) + body_3`
is equivalent to :code:`body_1 + (body_2 + body_3)`.
It is also commutative, up to some internal details which are usually not relevant.
However for more than two bodies, it is recommended to use instead the
:code:`join_bodies` method::

    all_bodies = body_1.join_bodies(body_2, body_3, body_4)

When two floating bodies with dofs are merged, the resulting body inherits from
the dofs of the individual bodies with the new name :code:`body_name__dof_name`.
For instance::

    body_1.add_translation_dof(name="Heave")
    body_2.add_translation_dof(name="Heave")
    both_bodies = body_1 + body_2
    assert 'body_1__Heave' in both_bodies.dofs
    assert 'body_2__Heave' in both_bodies.dofs


Clipping
--------

Meshes and bodies can be clipped with the :code:`clip` and :code:`clipped` methods.
As for the geometric transformations, the former is in-place whereas the second
returns a new object.
These methods take a :class:`~capytaine.meshes.geometry.Plane`
object as argument. The plane is defined by a point belonging to it and a normal
vector::

    xOy_Plane = Plane(point=(0, 0, 0), normal=(0, 0, 1))
    clipped_body = body.clipped(xOy_Plane)

Beware that the orientation of the normal vector of the :code:`Plane` will
determine which part of the mesh will be returned::

    higher_part = body.clipped(Plane(point=(0, 0, 0), normal=(0, 0, -1)))
    lower_part = body.clipped(Plane(point=(0, 0, 0), normal=(0, 0, 1)))
    # body = lower_part + higher_part

The method :code:`keep_immersed_part` will clip the body (by default in-place)
with respect to two horizontal planes at :math:`z=0` and :math:`z=-h`::

    clipped_body = body.keep_immersed_part(sea_bottom=-10, inplace=False)


Center of mass and rotation dofs
--------------------------------

The center of gravity of the body can be defined by assigning a vector of 3
elements to the :code:`center_of_mass` attribute::

    body.center_of_mass = np.array([0.0, -1.0, -1.0])

The center of mass is used in some hydrostatics computation.
It is not required for hydrodynamical coefficients, except for the definition of the rotation degrees of freedom.
When defining a rotation dof, the code looks for attributes called
:code:`rotation_center`, :code:`center_of_mass` or * :code:`geometric_center` (in that order),
and use them to define the rotation axis.
If none of them are define, the rotation is defined around the origin of the domain :math:`(0, 0, 0)`.


Defining an integration quadrature
----------------------------------

.. warning:: This feature is experimental.
             Only quadrilaterals panels are supported at the moment.

During the resolution of the BEM problem, the Green function has to be
integrated on the mesh. By default, the integration is approximated by taking
the value at the center of the panel and multiplying by its area. For a more
accurate intagration, an higher order quadrature can be defined.

This feature relies on the external package `quadpy` to compute the quadrature.
You can install it with::

    pip install quadpy

Then chose one of the `available quadratures
<https://github.com/nschloe/quadpy#quadrilateral>`_ and give it to the
:code:`compute_quadrature` method::

    from quadpy.quadrilateral import stroud_c2_7_2

    body.mesh.compute_quadrature(method=stroud_c2_7_2())

It will then be used automatically when needed.

.. warning:: Transformations of the mesh (merging, clipping, ...) may reset the quadrature.
             Compute it only on your final mesh.

