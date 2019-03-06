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
    import capytaine as cpt

    body = ...  # Set up the body and its dofs here.

    test_matrix = xr.Dataset(coords={
        'omega': np.linspace(0.1, 4, 40),
        'angle': [0, np.pi/2],
        'radiating_dof': list(body.dofs),
        'water_depth': [np.infty],
    })
    dataset = cpt.Nemoh().fill_dataset(test_matrix, [body])

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

If the coordinate :code:`theta` is added to the test matrix, the code will
compute the Kochin function for these values of :math:`\theta`.


The :class:`~capytaine.problems.LinearPotentialFlowProblem` class
-----------------------------------------------------------------

For a finer grain control of the problems to be solved, a list of :code:`LinearPotentialFlowProblem` should be defined.
It is defined as, e.g.::

    problem = cpt.RadiationProblem(body=my_body, omega=1.0, radiating_dof="Heave")
    other_problem = cpt.DiffractionProblem(body=my_body, omega=1.0, angle=pi/2)

Besides the body, all the parameters are optional.
The table below gives their definitions and their default values.

+-----------------------+------------------------------------------+------------------------+ 
| Parameter             | Description (unit)                       | Default value          |
+=======================+==========================================+========================+
| :code:`free_surface`  | Position of the free surface [#]_ (m)    | :math:`0` m            |
+-----------------------+------------------------------------------+------------------------+
| :code:`sea_bottom`    | Position of the sea bottom (m)           | :math:`-\infty` m      |
+-----------------------+------------------------------------------+------------------------+
| :code:`omega`         | Frequency :math:`\omega` (rad/s)         | :math:`1.0` rad/s      |
+-----------------------+------------------------------------------+------------------------+
| :code:`g`             | Acceleration of gravity :math:`g` (m/s²) | :math:`9.81` m/s²      |
+-----------------------+------------------------------------------+------------------------+
| :code:`rho`           | Water density (kg/m³)                    | :math:`1000` kg/m³     |
+-----------------------+------------------------------------------+------------------------+
| :code:`angle`         | Angle of incoming wave                   | :math:`0` rad          |
|                       | (only for diffraction)                   |                        |
+-----------------------+------------------------------------------+------------------------+
| :code:`radiating_dof` | Name of radiating dof                    | first one found        |
|                       | (only for radiation)                     |                        |
+-----------------------+------------------------------------------+------------------------+

.. [#] Only two positions are accepted for the free surface: :math:`z=0` and
       :math:`z= +\infty`. The latter corresponds to an object in a infinite
       domain filled with fluid and with no free surface.

The wave height is implicitely assumed to be :math:`1` m.
Since all computations are linear, any wave height or motion amplitude can be retrieved by multiplying the result by the desired value.

The following attributes are automatically computed for a given problem:

+------------------------------------+-------------------------------------------------+
| Parameter                          | Description (unit)                              |
+====================================+=================================================+
| :code:`depth`                      | Water depth :math:`h` (m)                       |
+------------------------------------+-------------------------------------------------+
| :code:`wavenumber`                 | Wave number :math:`k` (m¯¹)                     |
+------------------------------------+-------------------------------------------------+
| :code:`wavelength`                 | Wave length :math:`\lambda=\frac{2\pi}{k}` (m)  |
+------------------------------------+-------------------------------------------------+
| :code:`period`                     | Wave period :math:`T=\frac{2\pi}{\omega}` (s)   |
+------------------------------------+-------------------------------------------------+
| :code:`dimensionless_omega`        | :math:`\frac{2\omega^2 h}{g}` (ø)               |
+------------------------------------+-------------------------------------------------+
| :code:`dimensionless_wavenumber`   | :math:`k h` (ø)                                 |
+------------------------------------+-------------------------------------------------+

They can be retrieved as::

    problem.wavenumber
    problem.period
    # ...

Legacy Nemoh.cal parameters files
---------------------------------

The `legacy parameters files from Nemoh <https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-running-192930.kjsp>`_ can be read by a dedicated function::

    from capytaine.io.legacy import import_cal_file

    list_of_problems = import_cal_file("path/to/Nemoh.cal")

The function returns a list of :code:`LinearPotentialFlowProblems`.

.. warning:: This feature is experimental.
    Some of the settings in the files (such as the free surface computation or the Kochin function) are ignored for the moment.
    See the example :code:`Nemoh.cal` below.

.. literalinclude:: examples/Nemoh.cal

Command-line interface
----------------------

.. warning:: This feature is experimental.

Capytaine comes with a command-line command :code:`capytaine` which can be used as::

    $ capytaine path/to/directory/parameter.cal

The parameter file (in :code:`Nemoh.cal` format) passed as argument is read and legacy tecplot output file are written in the directory :code:`path/to/directory/results/`.

.. warning:: If results files already exist, they will be overwritten!

If no argument is provided to the command, the code looks for a file :code:`Nemoh.cal` in the current directory.

