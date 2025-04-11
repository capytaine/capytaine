========
Features
========

Capytaine is a Python library for the computation of wave loads on floating structures in frequency domain.
It can be used as a standard sea-keeping code for the standalone computation of hydrodynamical coefficients.
Its Python API also allows advanced users to integrate it in more complex workflows and to customize its behavior.

Capytaine is based on

* Nemoh_'s Fortran core routines for the computation of the Green function,
* meshmagick_ for the manipulation of meshes,
* xarray_ for the storage of hydrodynamical dataset,
* and various other tools from the `Python scientific ecosystem`_.

.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp
.. _meshmagick: https://github.com/LHEEA/meshmagick
.. _xarray: https://docs.xarray.dev
.. _`Python scientific ecosystem`: https://scipy.org/


Below some features of the latest version (version |version|) of Capytaine are listed.

Core features
-------------

* Computation of the added masses, radiation dampings, diffraction forces, Froude-Krylov forces and RAO.
* Single rigid body, multiple rigid bodies or bodies with any arbitrary degrees of freedom.
* Finite depth or deep water.
* Faster computation for plane-symmetric bodies.
* Approximate forward speed (for single rigid body only at the moment).
* Set up problems with either angular frequency, period, wavelength or angular wavenumber.
* Lid-based irregular frequencies removal.
* Computation of the pressure field, the velocity field, the free surface elevation and the Kochin function.
* Computation of the hydrostatic stiffness and rigid body inertia. Non-neutrally buoyant bodies are partially supported.
* 3D animations of the body motion and the free surface elevation.
* Various input mesh formats supported via meshmagick_.
* OpenMP parallelization.
* Output in NetCDF format (or legacy Nemoh Tecplot format).

Advanced features
-----------------

* A clean code with unit tests, meaningful variable names, comments, (almost) no Fortran 77, no global variables and no error message in French.
* Easy access to the internal of the solver through Python.
* Direct (potential formulation) or indirect (source formulation) boundary integral equation.
* Direct or iterative linear system solver.
* Single or double precision floating-point numbers.
* Several parameters to customize the evaluation of the Green function and its integration on the mesh.
* Possibility to plug-in other implementations of the Green function.

Experimental features
---------------------

* Hierarchical matrices with low-rank blocks for faster simulations.
* Toeplitz matrices for faster simulations for floating bodies with local symmetries or for regular arrays of identical floating bodies.

Roadmap
-------
See `on Github <https://github.com/orgs/capytaine/projects/1>`_.
