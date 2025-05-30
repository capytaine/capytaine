==========
Resolution
==========

Settings of the solver
----------------------
The :class:`~capytaine.bem.solver.BEMSolver` class takes three (keyword-only) arguments at the time of its initialization::

    from capytaine import BEMSolver
    solver = BEMSolver(green_function=..., engine=..., method=...)

Let us discuss in more details these parameters.

Green function
~~~~~~~~~~~~~~
A class used to evaluate the Green function, deriving from :class:`~capytaine.green_functions.abstract_green_function.AbstractGreenFunction`.
The following classes are available:

:class:`~capytaine.green_functions.delhommeau.Delhommeau` (Default)
   The method implemented in Nemoh (see [Del87]_ and [Del89]_).
   See the documentation for details on the available options.

   In Capytaine (and Nemoh), the integral of the wave term
   :math:`\mathcal{G}(r, z)` (and its derivative :math:`\frac{\partial
   \mathcal{G}}{\partial r}`) are approximated using surrogate models, which
   take the form of a tabulation of these function values for a grid of
   :math:`(r, z)`, precomputed at the initialization of the program. A
   third-order Lagrange polynomial interpolation is employed to obtain the
   values between the precomputed values.

   In version 1 of Capytaine (as in version 2 of Nemoh), the tabulation ranges
   of :math:`r` and :math:`z` are set as :math:`[0, 100]` with :math:`328`
   discretization values and :math:`[-16, 0]` with :math:`46` discretization
   values, respectively. In the new version, these can be user-defined with the
   following options::

        import capytaine as cpt

        # Legacy (versions 1.x)
        gf = cpt.Delhommeau(tabulation_nr=324, tabulation_rmax=100,
                            tabulation_nz=46, tabulation_zmin=-16,
                            tabulation_nb_integration_points=251,
                            tabulation_grid_shape="legacy",
                            finite_depth_method="legacy",
                            finite_depth_prony_decomposition_method="fortran",
                            gf_singularities="high_freq")

        # Default in Capytaine 2.1
        gf = cpt.Delhommeau(tabulation_nr=676, tabulation_rmax=100,
                            tabulation_nz=372, tabulation_zmin=-251,
                            tabulation_nb_integration_points=1001,
                            tabulation_grid_shape="scaled_nemoh3",
                            finite_depth_method="legacy",
                            finite_depth_prony_decomposition_method="fortran",
                            gf_singularities="high_freq")

        # Default in Capytaine 2.2
        gf = cpt.Delhommeau(tabulation_nr=676, tabulation_rmax=100,
                            tabulation_nz=372, tabulation_zmin=-251,
                            tabulation_nb_integration_points=1001,
                            tabulation_grid_shape="scaled_nemoh3",
                            finite_depth_method="legacy",
                            finite_depth_prony_decomposition_method="fortran",
                            gf_singularities="low_freq")

        # Default in Capytaine 2.3
        gf = cpt.Delhommeau(tabulation_nr=676, tabulation_rmax=100,
                            tabulation_nz=372, tabulation_zmin=-251,
                            tabulation_nb_integration_points=1001,
                            tabulation_grid_shape="scaled_nemoh3",
                            finite_depth_method="newer",
                            finite_depth_prony_decomposition_method="python",
                            gf_singularities="low_freq")

   In version 2.1, the default numbers of :math:`r` and :math:`z` values have
   been increased to :math:`676` and :math:`372`, respectively. While the range
   of :math:`r` is kept the same, the z range has been extended to
   :math:`[-251, 0]`. The option :code:`tabulation_grid_shape` is used to switched
   between the new distribution of points inspired by Nemoh version 3 or the
   :code:`"legacy"` approach. The :code:`tabulation_nb_integration_points`
   controls the accuracy of the precomputed tabulation points themselves.

   In version 2.2, the way singularities are extracted of the infinite depth
   Green function to be integrated has changed. The ``"low_freq"`` variant is
   expected to be more accurate at low frequency and near the free surface. The
   former variant is still available by setting the ``gf_singularities``
   parameter as in the above example.

   In version 2.3, some better variants to compute the finite depth Green
   function were introduced.

   The first time it is initialize with a given set of parameters, some tabulated
   data are precomputed and stored on disk.
   The default location is a os-dependant cache directory.
   The location at which the data is stored can be configured by passing
   ``tabulation_cache_dir`` to
   :class:`~capytaine.green_functions.delhommeau.Delhommeau` or by setting the
   environment variable ``CAPYTAINE_CACHE_DIR``.


:class:`~capytaine.green_functions.hams.LiangWuNoblesseGF`
   The infinite depth Green function from the following papers:

   [1] H. Wu, C. Zhang, Y. Zhu, W. Li, D. Wan, F. Noblesse,
       **A global approximation to the Green function for
       diffraction radiation of water waves**,
       Eur. J. Mech. B Fluids 65 (2017) 54-64.

   [2] H. Liang, H. Wu, F. Noblesse,
       **Validation of a global approximation for
       wave diffraction-radiation in deep water**,
       Appl. Ocean Res. 74 (2018) 80-86.

   Please cite them if you use this implementation.


:class:`~capytaine.green_functions.hams.FinGreen3D`
   The finite depth Green function from the following paper, as implemented in HAMS:

   Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun,
   Junliang Gao, Peiwen Cong, Guanghua He.
   **A reliable open-source package for performance evaluation of floating
   renewable energy systems in coastal and offshore regions**.
   Energy Conversion and Management, 174 (2018): 516-536.

   Please cite this paper if you use this implementation.


:class:`~capytaine.green_functions.hams.HAMS_GF`
   This class is just a thin wrapper around the two implementation above, using
   one or the other depending of the water depth.


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

   :code:`linear_solver` (Default: :code:`'lu_decomposition'`)
           This option is used to set the solver for linear systems that is used in the resolution of the BEM problem.
           Passing a string will make the code use one of the predefined solver. Three of them are available:
           :code:`'direct'` for a simple direct solver,
           :code:`'lu_decomposition'` for a faster direct solver with caching of the LU decomposition,
           or :code:`'gmres'` for an iterative solver.

           A direct solver is used by default (since version 1.4) because it is more robust and the computation time is more predictable.
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

Method
~~~~~~

The argument :code:`method` (default value: :code:`"indirect"`) controls
the approach employed to solve for the potential velocity solutions.
Two methods are implemented:

#. direct method (also known as "potential formulation", among other names)
   with :code:`method="direct"`,
#. indirect method (also known as "source formulation"), by default and with
   :code:`method="indirect"`.

The direct method appears to be slightly more accurate on some
test cases (especially when thin plates are involved) but is only implemented
for the computation of the forces on the floating body without forward speed.
Any other post-processing (e.g. free surface elevation) and forward speed
currently require the indirect method.

Since v2.3, the method is a parameter of :class:`~capytaine.bem.solver.BEMSolver`.
For backward compatibility, it can also be passed to
:meth:`~capytaine.bem.solver.BEMSolver.solve`,
:meth:`~capytaine.bem.solver.BEMSolver.solve_all` and
:meth:`~capytaine.bem.solver.BEMSolver.fill_dataset`, then overriding the
general setting of the solver.


Solving the problem
-------------------

Once the solver has been initialized, it can be used to solve problems with the
:meth:`~capytaine.bem.solver.BEMSolver.solve` method::

	result = solver.solve(problem, keep_details=False)

The optional argument :code:`keep_details` (default value: :code:`True`)
controls whether the source and potential distributions should be saved in the
result object. These data are necessary for some post-processing such as the
computation of the Kochin function or the reconstruction of the free surface
elevation. However, when only the force on the body is of interest, they can be
discarded to save space in memory.


A list of problems can be solved at once in an optimal order with::

	list_of_results = solver.solve_all(list_of_problems, keep_details=False)

where :meth:`~capytaine.bem.solver.BEMSolver.solve_all` accepts the same
optional keyword arguments as :meth:`~capytaine.bem.solver.BEMSolver.solve`.

When using :meth:`~capytaine.bem.solver.BEMSolver.solve_all`, a single problem
raising an error do not interrupt the full resolution. Instead, the error is
displayed in the log and the output result is replaced by a
:class:`~capytaine.bem.problems_and_results.FailedDiffractionResult` or a
:class:`~capytaine.bem.problems_and_results.FailedRadiationResult`.

Progress bar
------------

The methods :meth:`~capytaine.bem.solver.BEMSolver.solve_all` and
:meth:`~capytaine.bem.solver.BEMSolver.fill_dataset` display by default an
animated progress bar while solving.
This behavior can be turned off by giving the optional argument
``progress_bar=False`` to either method or by setting the environment variable
``CAPYTAINE_PROGRESS_BAR`` to ``False``.
This might be useful in testing environments and CI.

Timer
-----

The solver :class:`~capytaine.bem.solver.BEMSolver` keeps track of the time spent in some step of the resolution.
Results are stored in ``timer`` attribute and can also be accessed by :meth:`~capytaine.bem.solver.BEMSolver.timer_summary`.


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
