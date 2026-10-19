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

This example solves a BEM problem and display the computed pressure on the hull
on a 3D view of the mesh.

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

This example computes the free surface elevation as a post-processing step of
a diffraction problem.

.. literalinclude:: src/B3_free_surface_elevation.py
   :language: python


B4. Kochin function
-------------------

This example computes the Kochin function for a surging buoy and plot the
results.

.. literalinclude:: src/B4_kochin.py
   :language: python


B5. Plot velocity in domain
---------------------------

This example reconstruct the velocity field in the domain and display it
using matplotlib.

.. literalinclude:: src/B5_plot_velocity_in_domain.py
   :language: python


B6. Animated free surface elevation
-----------------------------------

This example solves a diffraction problem, it computes the free surface
elevation and shows it as a 3D animation.

.. literalinclude:: src/B6_animate_free_surface.py
   :language: python


B7. Animation of the RAO
------------------------

This script generates the animation of the RAO motion for a wave incoming in front of a ship,
such as the one used on the main page of this documentation.
This script requires the mesh of the ship :code:`boat_200.mar`. It can be
downloaded from: `<https://raw.githubusercontent.com/capytaine/capytaine/master/docs/examples/src/boat_200.mar>`_

.. literalinclude:: src/B7_boat_animation.py
   :language: python

B8. Pressure field for zero or infinite frequency
-------------------------------------------------

The computations for zero and infinite frequencies are tricky, because several magnitudes such as the pressure are zero or infinity, but they can be expressed as :math:`\omega` or :math:`\omega^2` times some finite value.
The example below shows an example of extracting this value.

.. literalinclude:: src/B8_pressure_infinite_frequency.py
   :language: python

B9. Mean drift force
--------------------

This example computes the mean drift force on a floating box with both the near field and the far field formulations, and compares them.
The near field mean drift force is still an experimental feature, whose interface may change in the future.

.. literalinclude:: src/B9_mean_drift_force.py
   :language: python

B10. Mean drift pressure field
------------------------------

This example extracts the pressure fields from the near field mean drift force computation, for a sea made of two waves with the same frequency:
the first-order pressure on the hull, the mean drift pressure on the hull, and the waterline term spread along a semi-elliptic profile around the mean waterline.
They are plotted on the mesh of the whole box, and the first-order pressure is animated over a wave period.
To save the animation as a video, uncomment the call to ``anim.save`` (requires `ffmpeg <https://ffmpeg.org/>`_).

.. literalinclude:: src/B10_second_order_pressure_field.py
   :language: python
