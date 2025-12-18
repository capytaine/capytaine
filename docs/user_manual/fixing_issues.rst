Common issues
=============

There are NaNs in my dataset output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Capytaine might leave a NaN value in the output dataset if the corresponding problem could not be solved or if it is ill-defined.

Normally, a message explaining the reason should be printed in the log.
You might need to turn on a more verbose logging to be sure to see it::

    cpt.set_logging('INFO')

In normal usage, the following cases cause a NaN:

* In low frequency and low depth ( :math:`k h < 0.1` ), then added mass, radiation damping and diffraction force are NaN . This is a limitation of the current implementation of the Green function used by default in Capytaine, and it might be improved in the future.

* In the limit cases :math:`\omega = 0` or :math:`\omega = \infty`:
    * the diffraction force and Froude-Krylov force are NaN, because they are ill-defined.
    * the added mass and radiation damping are always defined (not NaN) for :math:`\omega = \infty`
    * the added mass and radiation damping are defined for :math:`\omega = 0` in infinite depth, but they are NaN in finite depth, due to the same limitation of the implementation mentioned above.

* With forward speed, if the encounter frequency is zero, all outputs are NaN.
