==========================
Meshes and floating bodies
==========================

Importing a mesh
----------------

To create a new body using an existing mesh file, use the following syntax::

    from capytaine import FloatingBody

    body = FloatingBody.from_file('path/to/mesh.dat', file_format='nemoh')

The above example uses `Nemoh's mesh format`_.

.. _`Nemoh's mesh format`: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp

Thanks to Meshmagick, numerous other mesh format can be imported.
The file format can be given with the :code:`file_format` optional argument.
If no format is given, the code will try to infer it from the file extension::

    body = FloatingBody.from_file('path/to/mesh.msh')  # gmsh file

The formats currently supported by Meshmagick in reading are the following (from the documentation of Meshmagick):

+-----------+-----------------+----------------------+
| File      | Software        | Keywords             |
| extension |                 |                      |
+===========+=================+======================+
|   .mar    | NEMOH [#f1]_    | nemoh, mar           |
+-----------+-----------------+----------------------+
|   .nem    | NEMOH [#f1]_    | nemoh_mesh, nem      |
+-----------+-----------------+----------------------+
|   .gdf    | WAMIT [#f2]_    | wamit, gdf           |
+-----------+-----------------+----------------------+
|   .inp    | DIODORE [#f3]_  | diodore-inp, inp     |
+-----------+-----------------+----------------------+
|   .DAT    | DIODORE [#f3]_  | diodore-dat          |
+-----------+-----------------+----------------------+
|   .hst    | HYDROSTAR [#f4]_| hydrostar, hst       |
+-----------+-----------------+----------------------+
|   .nat    |    -            | natural, nat         |
+-----------+-----------------+----------------------+
|   .msh    | GMSH 2 [#f5]_   | gmsh, msh            |
+-----------+-----------------+----------------------+
|   .rad    | RADIOSS         | rad, radioss         |
+-----------+-----------------+----------------------+
|   .stl    |    -            | stl                  |
+-----------+-----------------+----------------------+
|   .vtu    | PARAVIEW [#f6]_ | vtu                  |
+-----------+-----------------+----------------------+
|   .vtp    | PARAVIEW [#f6]_ | vtp                  |
+-----------+-----------------+----------------------+
|   .vtk    | PARAVIEW [#f6]_ | paraview-legacy, vtk |
+-----------+-----------------+----------------------+
|   .tec    | TECPLOT [#f7]_  | tecplot, tec         |
+-----------+-----------------+----------------------+
|   .med    | SALOME [#f8]_   | med, salome          |
+-----------+-----------------+----------------------+

.. [#f1] NEMOH is an open source BEM Software for seakeeping developped at
         Ecole Centrale de Nantes (LHEEA)
.. [#f2] WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
.. [#f3] DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
.. [#f4] HYDROSTAR is a BEM Software for seakeeping developped by
         BUREAU VERITAS
.. [#f5] GMSH is an open source meshing software developped by C. Geuzaine
         and J.-F. Remacle. Version 4 of the file format is not supported at the
         moment.
.. [#f6] PARAVIEW is an open source visualization software developped by
         Kitware
.. [#f7] TECPLOT is a visualization software developped by Tecplot
.. [#f8] SALOME-MECA is an open source software for computational mechanics
         developped by EDF-R&D


Display
-------
Use the :code:`show` method to display the body in 3D using VTK::

    body.show()


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

