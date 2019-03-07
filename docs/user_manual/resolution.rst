==========
Resolution
==========

Settings of the solver
----------------------
Several parameters can be passed to the BEM solver at the
time of its initialization, such as::

    from capytaine import Nemoh
    solver = Nemoh(linear_solver="gmres", matrix_cache_size=0)

The most important of these parameters are described below.
See :class:`~capytaine.bem.nemoh.Nemoh` for a full list of accepted parameters.

:code:`hierarchical_matrices` (Default: :code:`True`)
	If :code:`True` and if the mesh of the body has a hierarchical structure,
	the solver will use this structure to build a Hierarchical Toeplitz matrix
	to save computation time. For advanced users, the solver parameters
	:code:`ACA_distance` and :code:`ACA_tol` can be use to set the precision of
	the Adaptive Cross Approximation.

:code:`linear_solver` (Default: :code:`'gmres'`)
	This option is used to set the solver for linear systems that is used in the resolution of the BEM problem.
	Passing a string will make the code use one of the predefined solver. Two of them are available:
	:code:`'direct'` for a direct solver using LU-decomposition or :code:`'gmres'` for an iterative solver.

	Alternatively, any function taking as arguments a matrix and a vector and returning a vector can be given to the solver::

		import numpy as np
		def my_linear_solver(A, b):
			"""A dumb solver for testing."""
			return np.linalg.inv(A) @ b
		my_bem_solver = Nemoh(linear_solver=my_linear_solver)

	This technique can be used for instance to apply a custom preconditioning to
	the iterative solver. It is recommended then to use
	:code:`hierarchical_matrices=False` to ensure that the custom function
	receives an usual numpy array.

:code:`matrix_cache_size` (Default: :code:`1`)
	The solver keeps in memory the last interaction matrices that has been computed.
	This setting controls the number of old matrices that are saved.
	Setting it to :code:`0` will reduce the RAM usage of the code but might
	increase the computation time.

:code:`cache_rankine_matrices` (Default: :code:`False`)
	If :code:`True`, the solve will cache separately the Rankine part of the
	influence matrix and the wave part. Indeed, since the former is independent
	of the wave frequency, it is not necessary to recompute it when studying the
	same mesh at different frequencies. However, the gain in computation time is
	low and this option might be conflicting with :code:`hierarchical_matrices`
	in some rare cases.


Solving the problem
-------------------

Once the solver has been initialized, it can be used to solve problems with the :meth:`~capytaine.bem.nemoh.Nemoh.solve` method::

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

The resolution should happen in parallel with OpenMP. The number of cores used
by OpenMP is controlled by the environment variables :code:`OMP_NUM_THREADS`
(for the computation of the Green function by capytaine itself) and
:code:`MKL_NUM_THREADS` (for the linear solver from Intel's MKL library
distributed with conda).

Accessing the influence matrices (for advanced users)
-----------------------------------------------------

To only compute the influence matrices, see the solver methods
:meth:`~capytaine.bem.nemoh.Nemoh.build_matrices`,
:meth:`~capytaine.bem.nemoh.Nemoh.build_matrices_rankine` and
:meth:`~capytaine.bem.nemoh.Nemoh.build_matrices_wave`

