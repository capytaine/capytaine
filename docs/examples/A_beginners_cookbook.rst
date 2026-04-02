===================
Beginner's cookbook
===================

This page contains several examples of Capytaine's main features, mostly about
the computation of added mass, radiation damping and excitation force.
The scripts can be downloaded individually as Python files from the |examples_folder|_.
See also the template in :doc:`../user_manual/quickstart`.

.. contents:: Content


A1. Single rigid body hydrodynamics
-----------------------------------

.. literalinclude:: src/A1_single_body_hydrodynamics.py
:language: python


A2. Simulation with several rigid bodies
----------------------------------------

.. literalinclude:: src/A2_multibody.py
   :language: python


A3. Finite depth flap
---------------------

This example illustrates the following features: finite depth, rotation around
an arbitrary axis (flap's hinge) and fixed obstacle without degrees of freedom.

.. literalinclude:: src/A3_finite_depth_flap.py
   :language: python


A4. Forward speed
-----------------

.. literalinclude:: src/A4_forward_speed_on_vertical_cylinder.py
   :language: python


A5. Benchmark plane symmetry
----------------------------

.. literalinclude:: src/A5_benchmark_plane_symmetries.py
   :language: python


A6. Benchmark axial symmetry
----------------------------

.. literalinclude:: src/A6_benchmark_axisymmetric_mesh.py
   :language: python


A7. Custom degree of freedom
----------------------------

This example defines arbitrary degrees of freedom for a sphere and solves a diffraction problem.

.. literalinclude:: src/A7_custom_dofs.py
   :language: python

The diffraction force on the "Heave" and "Bulge" dofs should be the same for both incoming wave directions.
The diffraction force on the "x-shear" dof is zero when the wave comes from the y direction.


A8. Convergence study
---------------------

This example runs a mesh convergence study for a submerged cylinder.

.. literalinclude:: src/A8_convergence_study.py
   :language: python


A9. Test irregular frequencies removal
--------------------------------------

This example compare the output without a lid and with a lid at different position.

.. literalinclude:: src/A9_test_irregular_frequency_removal.py
   :language: python


A10. Elastic beam
-----------------

.. literalinclude:: src/A10_elasticity_of_beam.py
   :language: python



A11. Parametric study: water depth
----------------------------------

This example runs the same simulation for several water depth and plot the results.

.. literalinclude:: src/A11_parametric_study_depth.py
   :language: python
