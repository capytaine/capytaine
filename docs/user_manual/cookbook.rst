Cookbook
========

This page contains several example of Capytaine's features.

.. contents:: Content


Added mass of a rigid body
--------------------------

This example generate the mesh of an horizontal cylinder, solves radiation problems for the six
rigid body degrees of freedom and then plot the computed added mass.

.. include:: examples/radiation_cylinder.py
   :code: python

Importing a mesh
----------------

Load a mesh from a file and display it with VTK.

.. include:: examples/import_mesh.py
   :code: python

The resulting body can be used to define a :code:`RadiationProblem` or a
:code:`DiffractionProblem`.

Symmetric body
--------------

This example loads a mesh from a file, keep only a part of it and define a symmetric body from this
half.

.. include:: examples/symmetric_body.py
   :code: python

Axisymmetric body
-----------------

This example generate an axisymmetric mesh from a profile function and solve radiation problems for
this floating body.

.. include:: examples/axisymmetric_buoy.py
   :code: python

Custom degree of freedom
------------------------

TODO

Simulation with several bodies
------------------------------

TODO

Animated free surface elevation
-------------------------------

This example solves a diffraction problem, it computes the free surface elevation and show it as a
3D animation.

.. include:: examples/animate_free_surface.py
   :code: python

Kochin function
---------------

TODO

Plot the influence matrix
-------------------------

This example plot the influence matrix for an horizontal cylinder.

.. include:: examples/plot_matrix.py
   :code: python

