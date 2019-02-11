======
Inputs
======

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
|   .msh    | GMSH [#f5]_     | gmsh, msh            |
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

.. rubric:: Footnotes

.. [#f1] NEMOH is an open source BEM Software for seakeeping developped at
         Ecole Centrale de Nantes (LHEEA)
.. [#f2] WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
.. [#f3] DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
.. [#f4] HYDROSTAR is a BEM Software for seakeeping developped by
         BUREAU VERITAS
.. [#f5] GMSH is an open source meshing software developped by C. Geuzaine
         and J.-F. Remacle
.. [#f6] PARAVIEW is an open source visualization software developped by
         Kitware
.. [#f7] TECPLOT is a visualization software developped by Tecplot
.. [#f8] SALOME-MECA is an open source software for computational mechanics
         developped by EDF-R&D

Transforming a mesh
-------------------

Several functions are available to transform existing meshes.

Most transformation methods exist in two versions: 

* one, named as a infinitive verb (`translate`, `rotate`, ...), is an in-place transformation;
* the other, named as a past participle (`translated`, `rotated`, ...), is the
  same transformation but returning a new object. 

In most cases, performance is not significant and the latter method should be
preferred since it makes code slightly easier to debug.

Below is a non-exhaustive list of available methods.
All of them can be applied to both meshes or floating bodies (although in the
latest version of the code, the implementation for floating bodies might not be
complete)::

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
    from capytaine.tools.geometry import Axis
    mesh.rotated(axis=Axis(vector=[1, 1, 1], point=[3, 4, 5]), angle=3.14/5)

    # Rotation such that the axis1 would become parallel to axis2
    mesh.rotated_to_align_axes(
        axis1=Axis(vector=[1, 1, 1]),
        axis2=Axis(vector=[3, 5, 7]),
        )

    # REFLECTIONS
    from capytaine.tools.geometry import Plane
    mesh.mirrored(Plane(normal=[1, 2, 1], point=[0, 4, 5]))


TODO: Document clipping functions.


Legacy Nemoh.cal parameters files
---------------------------------

The `legacy parameters files from Nemoh`_ can be read by a dedicated function::

    from capytaine.tools.import_export import import_cal_file

    list_of_problems = import_cal_file("path/to/Nemoh.cal")

.. _`legacy parameters files from Nemoh`: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-running-192930.kjsp

The function returns a list of :code:`LinearPotentialFlowProblems`.

.. warning:: This feature is experimental.
    Some of the settings in the files (such as the free surface computation or the Kochin function) are ignored for the moment.
    See the example :code:`Nemoh.cal` below.

.. literalinclude:: examples/Nemoh.cal

Command-line interface
----------------------

.. warning:: This feature is experimental.

Capytaine comes with a command-line command :code:`capytaine` which can be used as::

    $ capytaine path/to/directory/parameter.cal

The parameter file (in :code:`Nemoh.cal` format) passed as argument is read and legacy tecplot output file are written in the directory :code:`path/to/directory/results/`.

.. warning:: If results files already exist, they will be overwritten!

If no argument is provided to the command, the code looks for a file :code:`Nemoh.cal` in the current directory.


