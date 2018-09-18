====================
Setting up a problem
====================

Defining a test matrix with :code:`xarray`
------------------------------------------

The shortest way to set the parameters of the potential flow problems to be solved is to build a `xarray <http://xarray.pydata.org>`_ dataset.
Then the function :meth:`fill_dataset <capytaine.Nemoh.Nemoh.fill_dataset>` will automatically make all the simulations defined in the test matrix.

::

    import numpy as np
    import xarray as xr
    from capytaine import *

    body = ...  # Set up the body and its dofs here.

    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4, 40),
        'angle': [0, np.pi/2],
        'radiating_dof': list(body.dofs),
        'water_depth': [np.infty],
    })
    dataset = Nemoh().fill_dataset(test_matrix, [body])

It returns a dataset:

::

    >>> print(dataset)
    <xarray.Dataset>
    Dimensions:              (angle: 2, influenced_dof: 1, omega: 40, radiating_dof: 1)
    Coordinates:
      * omega                (omega) float64 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
      * radiating_dof        (radiating_dof) object 'Heave'
      * influenced_dof       (influenced_dof) object 'Heave'
      * angle                (angle) float64 0.0 1.571
    Data variables:
        added_mass           (omega, radiating_dof, influenced_dof) float64 2.14e+03 ...
        radiation_damping    (omega, radiating_dof, influenced_dof) float64 1.562e-06 ...
        Froude_Krylov_force  (omega, angle, influenced_dof) complex128 (-38.148887002448475-3.191891195797325e-16j) ...
        diffraction_force    (omega, angle, influenced_dof) complex128 (-21.354062211407268-1.5242955876404453e-07j) ...
    Attributes:
        g:            9.81
        rho:          1000.0
        body_name:    sphere_1
        water_depth:  inf


The :class:`~capytaine.problems.LinearPotentialFlowProblem` class
-----------------------------------------------------------------

.. note:: Work in progress...

For a finer grain control of the problems to be solved, a list of :code:`LinearPotentialFlowProblem` should be defined.

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
