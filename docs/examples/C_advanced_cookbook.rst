=================
Advanced cookbook
=================

This page contains several examples of Capytaine's advanced features, about the
fine control of the internal of the solver.
The scripts can be downloaded individually as Python files from the |examples_folder|_.

.. contents:: Content

C1. Parallel computation with Joblib
------------------------------------

TODO


C2. Direct or indirect boundary integral equation
-------------------------------------------------

TODO


C3. Single precision computations
---------------------------------

TODO


C4. Using an iterative solver
-----------------------------

TODO

C5. Plot the influence matrix
-----------------------------

This example plots the influence matrix for an horizontal cylinder.

.. literalinclude:: src/C5_plot_influence_matrix.py
   :language: python


C6. Toeplitz matrix of an axisymmetric buoy
-------------------------------------------

.. literalinclude:: src/C6_axisymmetric_buoy.py
   :language: python


C7. Hierarchical matrices with preconditionning
-----------------------------------------------

.. literalinclude:: src/C7_h_matrices_with_preconditionner.py
   :language: python



C8. Compare two implementations of the Green function
-----------------------------------------------------

This is an example of comparison of two implementations of the Green function.

.. literalinclude:: src/C8_compare_Green_functions.py
   :language: python


C9. Use a custom Green function
-------------------------------

This is an example of how to implement a custom Green function.

.. literalinclude:: src/C9_custom_Green_function.py
   :language: python


C10. Implementation of a GPU Custom Linear Solver
-------------------------------------------------

This is an example of how to implement a custom solver for linear systems that
leverages your computers's GPU, here with PyTorch.
Your mileage may vary as performance will vary wildly across different hardware.
Please also note that you are expected to already be familiar with GPU
frameworks and have a working setup to use this; Capytaine's developers can
offer no support for GPU issues unrelated to Capytaine.

.. literalinclude:: src/C10_custom_linear_solver_on_gpu.py
   :language: python
