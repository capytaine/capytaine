==========
Resolution
==========

Settings of the solver
----------------------
The :class:`~capytaine.bem.solver.BEMSolver` class takes two (keyword-only) arguments at the time of its initialization::

    from capytaine import BEMSolver
    solver = BEMSolver(green_function=..., engine=...)

Let us discuss in more details these two objects.

Green function
~~~~~~~~~~~~~~
A class used to evaluate the Green function, deriving from :class:`~capytaine.green_functions.abstract_green_function.AbstractGreenFunction`.
Two of them are available in the present version:

:class:`~capytaine.green_functions.delhommeau.Delhommeau` (Default)
   The method implemented in Nemoh (see [Del87]_ and [Del89]_).
   See the documentation for details on the available options.

:class:`~capytaine.green_functions.delhommeau.XieDelhommeau`
   A variant of the above, more accurate near the free surface (see [X18]_).
   Accepts the same options as :class:`Delhommeau <capytaine.green_functions.delhommeau.Delhommeau>`

Advanced users can write their own class to evaluate the Green function.
See the example in the :doc:`cookbook`.

Engine
~~~~~~
A class to build a interaction matrix, deriving from :class:`MatrixEngine <capytaine.bem.engines.MatrixEngine>`.
Two of them are available in the present version:

:class:`~capytaine.bem.engines.BasicMatrixEngine` (Default)
   A simple engine fairly similar to the one in Nemoh.
   It builds the full matrices with few optimizations.
   Only a reflection symmetry can be used to make the resolution faster.

   The object can be initialized with the following options:

   :code:`matrix_cache_size` (Default: :code:`1`)
           The solver keeps in memory the last interaction matrices that has been computed.
           This setting controls the number of old matrices that are saved.
           Setting it to :code:`0` will reduce the RAM usage of the code but might
           increase the computation time.

   :code:`linear_solver` (Default: :code:`'direct'`)
           This option is used to set the solver for linear systems that is used in the resolution of the BEM problem.
           Passing a string will make the code use one of the predefined solver. Two of them are available:
           :code:`'direct'` for a direct solver using LU-decomposition or :code:`'gmres'` for an iterative solver.

           The former is used by default (since version 1.4) because it is more robust and the computation time is more predictable.
           Advanced users might want to change the solver to :code:`gmres`, which is faster in many situations (and completely fails in other).

           Alternatively, any function taking as arguments a matrix and a vector and returning a vector can be given to the solver::

                   import numpy as np

                   def my_linear_solver(A, b):
                           """A dumb solver for testing."""
                           return np.linalg.inv(A) @ b

                   my_bem_solver = cpt.BEMSolver(
                      engine=BasicMatrixEngine(linear_solver=my_linear_solver)
                      )

           This option can be used for instance to apply a custom preconditioning to
           the iterative solver.

:class:`~capytaine.bem.engines.HierarchicalToeplitzMatrixEngine`
   Experimental engine using hierarchical structure in the mesh to build
   hierarchical influence matrices.

   The object can be initialized with the following options:

   :code:`matrix_cache_size` (Default: :code:`1`)
      Same as above.

   :code:`ACA_distance` and :code:`ACA_tol`
      Parameters of the Adaptive Cross Approximation (ACA) used to set the
      precision of the low-rank matrices.


Legacy interface
----------------

The class :class:`~capytaine.bem.solver.Nemoh` was the main solver class in
version 1.0 of Capytaine.
It is still available in the current version for backward compatibility.
It is now a subclass of :class:`~capytaine.bem.solver.BEMSolver` that always uses
:class:`~capytaine.green_functions.delhommeau.Delhommeau`'s Green function and
accept the same arguments as in version 1.0.

The use of :class:`~capytaine.bem.solver.BEMSolver` is recommended.

Solving the problem
-------------------

Once the solver has been initialized, it can be used to solve problems with the :meth:`~capytaine.bem.solver.BEMSolver.solve` method::

	result = solver.solve(problem, keep_details=False)

The optional argument :code:`keep_details` (default value: :code:`True`)
controls whether the source and potential distributions should be saved in the
result object. These data are necessary for some post-processing such as the
computation of the Kochin function or the reconstruction of the free surface
elevation. However, when only the force on the body is of interest, they can be
discarded to save space in memory.

A list of problems can be solved at once in an optimal order with::

	list_of_results = solver.solve_all(list_of_problems, keep_details=False)

Parallelization
---------------

Capytaine includes two kinds of parallelization.

+---------------------------+----------------+--------+
|                           | `joblib`       | OpenMP |
+---------------------------+----------------+--------+
| Single resolution         | ✗              | ✓      |
| (:code:`solve`)           |                |        |
+---------------------------+----------------+--------+
| Batch resolution          | ✓              | ✓      |
| (:code:`solve_all`        | (if installed) |        |
| and :code:`fill_dataset`) |                |        |
+---------------------------+----------------+--------+

Single problem with OpenMP
~~~~~~~~~~~~~~~~~~~~~~~~~~

When solving a single problem, matrix constructions and linear algebra
operations (using BLAS or MKL depending on your installation) can be
parallelized by OpenMP. This feature is installed and on by default. The number
of threads used can be controlled by the environment variable
:code:`OMP_NUM_THREADS`, as well as :code:`MKL_NUM_THREADS` (for the linear
algebra when using Intel's MKL library usually distributed with conda). Note
that the environment variable should be set *before* the start of the Python
interpreter. Alternatively, if you'd like to change dynamically the number of
threads, it can be done with the `threadpoolctl library
<https://github.com/joblib/threadpoolctl>`_ (see also :issue:`47`).

Batch resolution with joblib
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When solving several independent problems, they can be solved in parallel. This
feature (new in version 1.4) requires the optional dependency `joblib
<https://github.com/joblib/joblib>`_ to be installed. The methods
:meth:`~capytaine.bem.solver.BEMSolver.solve_all` and
:meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` take an optional
keyword-argument :code:`n_jobs` which control the number of jobs to run in
parallel during the batch resolution.
Since `joblib` may disturb user feedback (logging and error
reporting), it is currently disabled by default.

When :code:`n_jobs=1` (the default) or `joblib` is not installed, no parallel
batch resolution happens (although OpenMP parallelization might still be
enabled).

When :code:`n_jobs=-1`, all CPU cores are used (and `joblib` should
automatically disable the OpenMP parallelization.)

The two parallelization layers (OpenMP and `joblib`) have different usage. If
you have a relatively small mesh but study a large number of sea states, you
should use the `joblib` parallelization. On the other hand, if your mesh is
large or your available RAM is low, it might be beneficial to turn off the
`joblib` parallelization and use only the OpenMP one.
