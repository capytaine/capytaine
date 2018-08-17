====================
Setting up a problem
====================

.. note:: Work in progress...

The :code:`LinearPotentialFlowProblem` class
--------------------------------------------

Main parameters
~~~~~~~~~~~~~~~

+---------------+------------------------------------------+------------------------+ 
| Parameter     | Description (unit)                       | Default value          |
+===============+==========================================+========================+
| free_surface  | Position of the free surface (m)         | :math:`0` m            |
+---------------+------------------------------------------+------------------------+
| sea_bottom    | Position of the sea bottom (m)           | :math:`-\infty` m      |
+---------------+------------------------------------------+------------------------+
| omega         | Frequency :math:`\omega` (rad/s)         | :math:`1.0` rad/s      |
+---------------+------------------------------------------+------------------------+
| g             | Acceleration of gravity :math:`g` (m/s²) | :math:`9.81` m/s²      |
+---------------+------------------------------------------+------------------------+
| rho           | Water density (kg/m³)                    | :math:`1000` kg/m³     |
+---------------+------------------------------------------+------------------------+
| angle         | Angle of incoming wave                   | :math:`0` rad          |
|               | (only for diffraction)                   |                        |
+---------------+------------------------------------------+------------------------+
| radiating_dof | Name of radiating dof                    | first one found        |
|               | (only for radiation)                     |                        |
+---------------+------------------------------------------+------------------------+

The wave height is implicitely assumed to be :math:`1` m.
Since all computations are linear, any wave height or motion amplitude can be retrieved by multiplying the result by the desired value.

Derived parameters
~~~~~~~~~~~~~~~~~~

+----------------------------+-------------------------------------------------+
| Parameter                  | Description (unit)                              |
+============================+=================================================+
| depth                      | Water depth :math:`h` (m)                       |
+----------------------------+-------------------------------------------------+
| wavenumber                 | Wave number :math:`k` (m¯¹)                     |
+----------------------------+-------------------------------------------------+
| wavelength                 | Wave length :math:`\lambda=\frac{2\pi}{k}` (m)  |
+----------------------------+-------------------------------------------------+
| period                     | Wave period :math:`T=\frac{2\pi}{\omega}` (s)   |
+----------------------------+-------------------------------------------------+
| dimensionless_omega        | :math:`\frac{2\omega^2 h}{g}` (ø)               |
+----------------------------+-------------------------------------------------+
| dimensionless_wavenumber   | :math:`k h` (ø)                                 |
+----------------------------+-------------------------------------------------+
