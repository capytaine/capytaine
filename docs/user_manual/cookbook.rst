========
Cookbook
========

This page contains several examples of Capytaine's features.

.. contents:: Content

Base examples
=============

Importing a mesh
----------------

This example loads a mesh from a file and displays it with VTK.

.. literalinclude:: examples/import_mesh.py
   :language: python

The resulting body can be used to define a :code:`RadiationProblem` or a
:code:`DiffractionProblem`.

Custom degree of freedom
------------------------

This example defines arbitrary degrees of freedom for a sphere and solves a diffraction problem.

.. literalinclude:: examples/custom_dofs.py
   :language: python

The diffraction force on the "Heave" and "Bulge" dofs should be the same for both incoming wave directions.
The diffraction force on the "x-shear" dof is zero when the wave comes from the y direction.

Added mass of a rigid body
--------------------------

This example generates the mesh of an horizontal cylinder, solves radiation problems for the six
rigid-body degrees of freedom and then plots the added mass.

.. literalinclude:: examples/radiation_cylinder.py
   :language: python

Influence of the water depth
----------------------------

This example runs the same simulation for several water depth and plot the results.

.. literalinclude:: examples/finite_depth_cylinder.py
   :language: python

Simulation with several bodies
------------------------------

TODO


Intermediate examples
=====================

Symmetric body
--------------

This example loads a mesh from a file, keeps only a part of it and defines a symmetric body from this
half.

.. literalinclude:: examples/symmetric_body.py
   :language: python

Axisymmetric body
-----------------

This example generates an axisymmetric mesh from a profile function and solves radiation problems for
this floating body.

.. literalinclude:: examples/axisymmetric_buoy.py
   :language: python

Animated free surface elevation
-------------------------------

This example solves a diffraction problem, it computes the free surface elevation and shows it as a
3D animation.

.. literalinclude:: examples/animate_free_surface.py
   :language: python

Kochin function
---------------

TODO

Advanced examples
=================

Plot the influence matrix
-------------------------

This example plots the influence matrix for an horizontal cylinder.

.. literalinclude:: examples/plot_matrix.py
   :language: python

Parallel computations
---------------------

TODO

