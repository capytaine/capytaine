========
Features
========

Capytaine is a Python library for the computation of linear wave loads on floating structures in frequency domain.
It can be used as a standard sea-keeping code for the standalone computation of hydrodynamical coefficients.
Its Python user interface also allows advanced users to integrate it in more complex workflows and to customize its behavior.

Capytaine started as a Python interface to Nemoh_'s  core Fortran routines, interfacing them with the `Python scientific ecosystem`_.
Since then, it has evolved independently of Nemoh, and offers several other backends and a few features that are not found in Nemoh (which also received new features in the mean time with version 3).

.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp
.. _`Python scientific ecosystem`: https://scipy.org/

Below some features of the latest version (version |version|) of Capytaine are listed.

Core features
-------------

* Computation of the added masses, radiation dampings, diffraction forces, Froude-Krylov forces and RAO.
* Single rigid body, multiple rigid bodies or bodies with any generalized degrees of freedom.
* Finite depth or deep water.
* Approximate forward speed (for single rigid body only at the moment).
* Set up problems with either angular frequency, period, wavelength or angular wavenumber.
* Lid-based irregular frequencies removal.
* Post-processing computation of the pressure field, the velocity field, the free surface elevation and the Kochin function.
* Computation of the hydrostatic stiffness and rigid body inertia. Non-neutrally buoyant bodies are partially supported.
* Faster computation for plane-symmetric bodies.
* OpenMP threads parallelization, as well as optional processes parallelisation with joblib_.
* Output in NetCDF format through xarray_.
* 3D animations of the body motion and the free surface elevation.

.. _xarray: https://docs.xarray.dev
.. _joblib: https://github.com/joblib/joblib

Support tools
-------------

Besides its core features, Capytaine also offers support tools such as basic meshing tools to clip the mesh at the free surface or generate a lid.

For this kind of annex tools, the philosophy of Capytaine is to:
* Have built-in the bare minimal implementation to test the corresponding feature.
* Make it easy to interface with external software or Python libraries offering more advanced implementation for fine tuning. In this case, make it easy to interact with a meshing software to clip the mesh and generate a lid mesh.


Advanced features
-----------------

* A clean code with unit tests, meaningful variable names, comments, (almost) no Fortran 77, no global variables and no error message in French.
* Easy access to the internal of the solver through Python.
* Direct (potential formulation) or indirect (source formulation) boundary integral equation.
* Direct or iterative linear system solver.
* Single or double precision floating-point numbers.
* Several parameters to customize the evaluation of the Green function and its integration on the mesh.
* HAMS' Green function is available as a backend that can be used instead of Nemoh's one.
* Possibility to plug-in other implementations of the Green function.

Roadmap
-------
See `on Github <https://github.com/orgs/capytaine/projects/1>`_.
