============================================
Changing the performance/accuracy compromise
============================================

Identifying accuracy bottleneck and providing options to work around them is a current field of research around Capytaine.
See for instance `this paper <https://hal.science/hal-04822282>`_ and `the associated presentation<https://hal.science/hal-04822330/>`_ for some theoretical insights.


How to detect numerical errors
==============================

* Symmetry of the added mass and radiation damping matrices

* Mesh convergence test.

* Test on easy cases: 0 and infinite frequency problems in infinite depth are easy and should present much less numerical inaccuracies.


Parameters that can be changed
==============================

The default parameters are set to give a good compromise between performance and accuracy.
In this section, we list some of the parameter that users can change to try to increase the accuracy of the solver or to test the sensibility of their results.

No cost changes
~~~~~~~~~~~~~~~

* **Change the solver to direct solver.**
  The default indirect solver allows more post-processing and features, but the direct solver is often more accurate.
  The direct solver can be used by passing the ``method="direct"`` argument to `:meth:`~capytaine.bem.solver.BEMSolver.solve` or `:meth:`~capytaine.bem.solver.BEMSolver.solve_all` or `:meth:`~capytaine.bem.solver.BEMSolver.fill_dataset`::

    solver = cpt.BEMSolver()
    results = solver.solve_all(problems, method="direct")

  Computation time is the same for both solver.


* **Use the high-frequency asymptotics Green function.**
  The high-frequency asymptotics can be better captured by pssing the ``gf_singularities="high_freq"`` to the Green function::

    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="high_freq"))

  All other frequencies are usually more accurately computed with the default ``"low_freq"`` version.
  Computation time should be the same.


Slightly slower
~~~~~~~~~~~~~~~

* **Thiner tabulation.**

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


* **Mesh refinement.**
  Especially near the free surface and near the edges and corners of the mesh.
  There is no tool to do that directly in Capytaine at the moment.
  Todo: add reference on cosine spacing.


Probably not worth the computation cost, but you never know
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Disable tabulation of the Green function.**

* **Higher order quadrature scheme.**
