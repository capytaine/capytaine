========
Features
========

Capytaine is based on

* Nemoh_'s Fortran core routines for the computation of the Green function,
* meshmagick_ for the manipulation of meshes,
* and various tools from the `Python scientific ecosystem`_.

.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp
.. _meshmagick: https://github.com/LHEEA/meshmagick
.. _`Python scientific ecosystem`: https://scipy.org/

For users, it can be seen as a Python 3.6 (or higher) interface to Nemoh.
Since most of the code has been rewritten, it can be used by developpers as a
more concise and better documented version of Nemoh.

Below some features of the latest version (version |version|) of Capytaine are listed.

Main features
-------------

* Computation of the added masses, radiation dampings, diffraction forces and Froude-Krylov forces for rigid bodies or for bodies with any arbitrary degrees of freedom.
* Python object oriented API.
* Choice among several implementations of the Green function.
* OpenMP parallelization.
* Iterative linear system solver.
* Faster than Nemoh 2.0 on a single thread [#]_.
* Various input mesh formats supported via meshmagick_.
* Output in legacy Nemoh Tecplot format or NetCDF format.
* Computation of the free surface elevation and the Kochin function.
* 3D animations of the body motion and the free surface elevation.
* Binaries distributed via `conda <https://www.anaconda.com/download/>`_.
* A cleaner code with unit tests, meaningful variable names, more comments, no Fortran 77, no global variables and no error message in French.
* Double precision floating-point numbers by default. (Can be recompiled in single precision.)

Experimental features
---------------------

* Input via legacy ``Nemoh.cal`` files.
* Computation of the Response Amplitude Operators (RAO).
* Hierarchical matrices with low-rank blocks for faster simulations.
* Toeplitz matrices for faster simulations for floating bodies with local symmetries or for regular arrays of identical floating bodies.
* Higher order quadratures for more precise integration of the Green function.


.. rubric:: Footnotes

.. [#] However, the RAM usage of Capytaine is usually higher than the one of Nemoh 2.0.
   If your problem is memory bounded, simply switching to Capytaine won't help.
   But some advanced features such as the hierarchical matrices can help you reduce the cost of the problem.
