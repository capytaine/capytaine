===================
Beginner's cookbook
===================

This page contains several examples of Capytaine's main features, mostly about
the computation of added mass, radiation damping and excitation force.
The scripts can be downloaded individually as Python files from the |examples_folder|_.
See also the template in :doc:`../user_manual/quickstart`.

.. contents:: Content


A1. Added mass of a single rigid body
-------------------------------------

This example generates the mesh of an horizontal cylinder, solves radiation problems for the six
rigid-body degrees of freedom and then plots the added mass.

.. literalinclude:: src/A1_radiation_cylinder.py
   :language: python


A2. Simulation with several bodies
----------------------------------

.. literalinclude:: src/A2_multibody.py
   :language: python


A3. Finite depth flap
---------------------

This example illustrates the following features: finite depth, rotation around
an arbitrary axis (flap's hinge) and fixed obstacle without degrees of freedom.

.. literalinclude:: src/A3_finite_depth_flap.py
   :language: python


A4. Custom degree of freedom
----------------------------

This example defines arbitrary degrees of freedom for a sphere and solves a diffraction problem.

.. literalinclude:: src/A4_custom_dofs.py
   :language: python

The diffraction force on the "Heave" and "Bulge" dofs should be the same for both incoming wave directions.
The diffraction force on the "x-shear" dof is zero when the wave comes from the y direction.


A5. Convergence study
---------------------

This example runs a mesh convergence study for a submerged cylinder.

.. literalinclude:: src/A5_convergence_study.py
   :language: python


A5. Forward speed
-----------------

TODO


A6. Irregular frequencies removal
---------------------------------

.. literalinclude:: src/A6_irregular_frequency_removal.py
   :language: python


A7. Elastic beam
----------------

.. literalinclude:: src/A7_elasticity_of_beam.py
   :language: python


A8. Export dataset
------------------

.. literalinclude:: src/A8_export_dataset.py
   :language: python


A9. Paremetric study with several water depth
---------------------------------------------

This example runs the same simulation for several water depth and plot the results.

.. literalinclude:: src/A9_parametric_study_depth.py
   :language: python
