========
Cookbook
========

This page contains several examples of Capytaine's features.
The scripts can be downloaded individually as Python files from the |examples_folder|_.

.. contents:: Content

Base examples
=============

Added mass of a rigid body
--------------------------

This example generates the mesh of an horizontal cylinder, solves radiation problems for the six
rigid-body degrees of freedom and then plots the added mass.

.. literalinclude:: examples/radiation_cylinder.py
   :language: python

Custom degree of freedom
------------------------

This example defines arbitrary degrees of freedom for a sphere and solves a diffraction problem.

.. literalinclude:: examples/custom_dofs.py
   :language: python

The diffraction force on the "Heave" and "Bulge" dofs should be the same for both incoming wave directions.
The diffraction force on the "x-shear" dof is zero when the wave comes from the y direction.

Influence of the water depth
----------------------------

This example runs the same simulation for several water depth and plot the results.

.. literalinclude:: examples/finite_depth_cylinder.py
   :language: python


Irregular frequencies removal
-----------------------------

.. literalinclude:: examples/irregular_frequency_removal.py
   :language: python


Simulation with several bodies
------------------------------

.. literalinclude:: examples/multibody.py
   :language: python

Plot pressure on hull
---------------------

.. literalinclude:: examples/pressure_on_hull.py
   :language: python

Free surface elevation
----------------------

.. literalinclude:: examples/free_surface_elevation.py
   :language: python

Comparison with hydrostatics from Meshmagick
--------------------------------------------

Hydrostatic and inertia properties can be computed independently via Capytaine or Meshmagick.
This script compares them both with analytical expression for a simple geometric object.

.. literalinclude:: examples/hydrostatics.py
   :language: python

Intermediate examples
=====================

Animated free surface elevation
-------------------------------

This example solves a diffraction problem, it computes the free surface elevation and shows it as a
3D animation.

.. literalinclude:: examples/animate_free_surface.py
   :language: python

Animation of the RAO
--------------------

This script generates the animation of the RAO motion for a wave incoming in front of a ship,
such as the one used on the main page of this documentation.
This script requires the mesh of the ship :code:`boat_200.mar`. It can be
downloaded from: `<https://raw.githubusercontent.com/capytaine/capytaine/master/docs/user_manual/examples/boat_200.mar>`_

.. literalinclude:: examples/boat_animation.py
   :language: python

Kochin function
---------------

This example computes the Kochin function for a surging buoy and plot the results.

.. literalinclude:: examples/kochin.py
   :language: python

Haskind's relation
------------------

This example computes the excitation force from the radiation potential
using Haskind's relation. The result is compared with the one obtained by
direct integration of the potentials from incident waves and from the
diffraction problem.

.. literalinclude:: examples/haskind.py
    :language: python


Axisymmetric body
-----------------

This example generates an axisymmetric mesh from a profile function and solves radiation problems for
this floating body.

.. literalinclude:: examples/axisymmetric_buoy.py
   :language: python

Advanced examples
=================

Convergence study
-----------------

This example runs a mesh convergence study for a submerged cylinder.

.. literalinclude:: examples/convergence_study.py
   :language: python

Plot the influence matrix
-------------------------

This example plots the influence matrix for an horizontal cylinder.

.. literalinclude:: examples/plot_influence_matrix.py
   :language: python

Compare two implementations of the Green function
-------------------------------------------------

This is an example of comparison of two implementations of the Green function.

.. literalinclude:: examples/compare_Green_functions.py
   :language: python

Use a custom Green function
---------------------------

This is an example of how to implement a custom Green function.

.. literalinclude:: examples/custom_Green_function.py
   :language: python
