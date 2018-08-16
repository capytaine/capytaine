======
Inputs
======

Importing a mesh
----------------

To create a new body using an existing mesh file, use the following syntax::

    from capytaine import FloatingBody

    body = FloatingBody.from_file('path/to/mesh.mar')

The default format is `Nemoh's mesh format`_.

.. _`Nemoh's mesh format`: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-mesh-192932.kjsp

Thanks to Meshmagick, numerous other mesh format can be imported.
If your mesh file is in another format, you have to specify it by passing a keyword as an argument::

    body = FloatingBody.from_file('path/to/mesh.msh', file_format='gmsh')

The formats currently supported by Meshmagick are the following.

+-----------+------------+-----------------+----------------------+
| File      | R: Reading | Software        | Keywords             |
| extension | W: writing |                 |                      |
+===========+============+=================+======================+
|   .mar    |    R/W     | NEMOH [#f1]_    | nemoh, mar           |
+-----------+------------+-----------------+----------------------+
|   .nem    |    R       | NEMOH [#f1]_    | nemoh_mesh, nem      |
+-----------+------------+-----------------+----------------------+
|   .gdf    |    R/W     | WAMIT [#f2]_    | wamit, gdf           |
+-----------+------------+-----------------+----------------------+
|   .inp    |    R       | DIODORE [#f3]_  | diodore-inp, inp     |
+-----------+------------+-----------------+----------------------+
|   .DAT    |    W       | DIODORE [#f3]_  | diodore-dat          |
+-----------+------------+-----------------+----------------------+
|   .hst    |    R/W     | HYDROSTAR [#f4]_| hydrostar, hst       |
+-----------+------------+-----------------+----------------------+
|   .nat    |    R/W     |    -            | natural, nat         |
+-----------+------------+-----------------+----------------------+
|   .msh    |    R       | GMSH [#f5]_     | gmsh, msh            |
+-----------+------------+-----------------+----------------------+
|   .rad    |    R       | RADIOSS         | rad, radioss         |
+-----------+------------+-----------------+----------------------+
|   .stl    |    R/W     |    -            | stl                  |
+-----------+------------+-----------------+----------------------+
|   .vtu    |    R/W     | PARAVIEW [#f6]_ | vtu                  |
+-----------+------------+-----------------+----------------------+
|   .vtp    |    R/W     | PARAVIEW [#f6]_ | vtp                  |
+-----------+------------+-----------------+----------------------+
|   .vtk    |    R/W     | PARAVIEW [#f6]_ | paraview-legacy, vtk |
+-----------+------------+-----------------+----------------------+
|   .tec    |    R/W     | TECPLOT [#f7]_  | tecplot, tec         |
+-----------+------------+-----------------+----------------------+
|   .med    |    R       | SALOME [#f8]_   | med, salome          |
+-----------+------------+-----------------+----------------------+

.. rubric:: Footnotes

.. [#f1] NEMOH is an open source BEM Software for seakeeping developped at
         Ecole Centrale de Nantes (LHEEA)
.. [#f2] WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
.. [#f3] DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
.. [#f4] HYDROSTAR is a BEM Software for seakeeping developped by
         BUREAU VERITAS
.. [#f5] GMSH is an open source meshing software developped by C. Geuzaine
         and J.-_faces. Remacle
.. [#f6] PARAVIEW is an open source visualization software developped by
         Kitware
.. [#f7] TECPLOT is a visualization software developped by Tecplot
.. [#f8] SALOME-MECA is an open source software for computational mechanics
         developped by EDF-R&D


Legacy Nemoh.cal parameters files
---------------------------------

The legacy parameters files from Nemoh can be read by a dedicated function::

    from capytaine.tools.import_export import import_cal_file

    list_of_problems = import_cal_file("path/to/Nemoh.cal")

The function returns a list of :code:`LinearPotentialFlowProblems`.

.. warning:: This feature is experimental.
    Some of the settings in the files (such as the free surface computation or the Kochin function) are ignored for the moment.

