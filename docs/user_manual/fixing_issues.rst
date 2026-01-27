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


My results does not match the ones of this other code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Difference in conventions between Capytaine and other solvers are listed in :doc:`conventions`.


Besides, many other codes such as WAMIT and HAMS are using another formulation of the Boundary Integral Equations, which is implemented in Capytaine but not used by default.
You can turn it on by using the following in your script::

    solver = cpt.BEMSolver(method="direct")

The default ``"indirect"`` solver allows more post-processing and features, but the ``"direct"`` solver is often more accurate.
Computation time is the same for both solver.


If you believe the approximation method used by Capytaine to compute the Green function to be the cause of a discrepancy, you can refine it using the following::

  ::

    finer_gf = cpt.Delhommeau(
        tabulation_nr=4000, tabulation_rmax=100,
        tabulation_nz=2000, tabulation_zmin=-251,
        tabulation_nb_integration_points=5_000,
        tabulation_grid_shape="scaled_nemoh3",
    )
    solver = cpt.BEMSolver(green_function=finer_gf)

Increasing the tabulation is expensive the first time you call it as a new
tabulation needs to be constructed, but subsequent resolutions should use the
same tabulation and be as cheap as before.


The added mass and radiation damping matrices are not symmetric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Indeed, the added mass and radiation damping are theoretically expected to be symmetric.
They are usually not in Capytaine results due to numerical accuracy: solving
the same problem in two different ways leads to slightly different
approximation errors.
If the mismatch is big, this is the sign that the numerical error is high.
The asymmetry should disappear asymptotically for fine meshes.
