=================
Advanced cookbook
=================

This page contains several examples of Capytaine's advanced features, about the
fine control of the internal of the solver.
The scripts can be downloaded individually as Python files from the |examples_folder|_.

.. contents:: Content

Parallel computation with Joblib
--------------------------------

TODO


Direct or indirect boundary integral equation
---------------------------------------------

TODO


Iterative solver with preconditionning
--------------------------------------

.. literalinclude:: src/C_preconditioner.py
   :language: python


Plot the influence matrix
-------------------------

This example plots the influence matrix for an horizontal cylinder.

.. literalinclude:: src/C_plot_influence_matrix.py
   :language: python


Toeplitz matrix of an axisymmetric buoy
---------------------------------------

.. literalinclude:: src/C_axisymmetric_buoy.py
   :language: python


Single precision computations
-----------------------------

TODO


Compare two implementations of the Green function
-------------------------------------------------

This is an example of comparison of two implementations of the Green function.

.. literalinclude:: src/C_compare_Green_functions.py
   :language: python


Use a custom Green function
---------------------------

This is an example of how to implement a custom Green function.

.. literalinclude:: src/C_custom_Green_function.py
   :language: python
