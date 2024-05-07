=====================
Intermediate cookbook
=====================

This page contains several examples of Capytaine's features, mostly focused on
post-processing features such as reconstruction of the free surface elevation
or computation of the Kochin functions.
The scripts can be downloaded individually as Python files from the |examples_folder|_.

.. contents:: Content


B1. Plot pressure on hull
-------------------------

.. literalinclude:: src/B1_pressure_on_hull.py
   :language: python


B2. Haskind's relation
----------------------

This example computes the excitation force from the radiation potential
using Haskind's relation. The result is compared with the one obtained by
direct integration of the potentials from incident waves and from the
diffraction problem.

.. literalinclude:: src/B2_haskind.py
    :language: python


B3. Free surface elevation
--------------------------

.. literalinclude:: src/B3_free_surface_elevation.py
   :language: python


B4. Kochin function
-------------------

This example computes the Kochin function for a surging buoy and plot the results.

.. literalinclude:: src/B4_kochin.py
   :language: python


B5. Plot velocity in domain
---------------------------


.. literalinclude:: src/B5_plot_velocity_in_domain.py
   :language: python


B6. Animated free surface elevation
-----------------------------------

This example solves a diffraction problem, it computes the free surface elevation and shows it as a
3D animation.

.. literalinclude:: src/B6_animate_free_surface.py
   :language: python


B7. Animation of the RAO
------------------------

This script generates the animation of the RAO motion for a wave incoming in front of a ship,
such as the one used on the main page of this documentation.
This script requires the mesh of the ship :code:`boat_200.mar`. It can be
downloaded from: `<https://raw.githubusercontent.com/capytaine/capytaine/master/docs/user_manual/src/boat_200.mar>`_

.. literalinclude:: src/B7_boat_animation.py
   :language: python
