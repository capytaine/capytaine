Cookbook
========

This page contains several example of Capytaine's features.

.. contents:: Content

Custom degree of freedom
------------------------

This example defines arbitrary degrees of freedom for a sphere and solves a diffraction problem.

.. include:: examples/custom_dofs.py
   :code: python

The force on the "Heave" and "Bulge" dofs should be the same for both incoming wave directions.
The force on the "x-shear" dof is zero when the wave comes from the y direction.

Added mass of a rigid body
--------------------------

This example generates the mesh of an horizontal cylinder, solves radiation problems for the six
rigid body degrees of freedom and then plots the computed added mass.

.. include:: examples/radiation_cylinder.py
   :code: python

Importing a mesh
----------------

This example loads a mesh from a file and displays it with VTK.

.. include:: examples/import_mesh.py
   :code: python

The resulting body can be used to define a :code:`RadiationProblem` or a
:code:`DiffractionProblem`.

Symmetric body
--------------

This example loads a mesh from a file, keeps only a part of it and defines a symmetric body from this
half.

.. include:: examples/symmetric_body.py
   :code: python

Axisymmetric body
-----------------

This example generates an axisymmetric mesh from a profile function and solves radiation problems for
this floating body.

.. include:: examples/axisymmetric_buoy.py
   :code: python

Simulation with several bodies
------------------------------

TODO

Animated free surface elevation
-------------------------------

This example solves a diffraction problem, it computes the free surface elevation and shows it as a
3D animation.

.. include:: examples/animate_free_surface.py
   :code: python

Kochin function
---------------

TODO

Plot the influence matrix
-------------------------

This example plots the influence matrix for an horizontal cylinder.

.. include:: examples/plot_matrix.py
   :code: python

