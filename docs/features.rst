========
Features
========

Capytaine is based on

* Nemoh_ Fortran core routines for the computation of the Green function,
* meshmagick_ for the manipulation of meshes,
* and various tools from the `Python scientific ecosystem`_.

.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp
.. _meshmagick: https://github.com/LHEEA/meshmagick
.. _`Python scientific ecosystem`: https://scipy.org/

For users, it can be seen as a Python [#]_ interface to Nemoh.
Since most of the code has been rewritten, it can be used by developpers as a
more concise and better documented version of Nemoh.

Below some features of the latest version (|version|) of Capytaine are listed.

Main features
-------------

* Computation of the **added masses, radiation dampings, diffraction forces and Froude-Krylov forces** for rigid bodies or for bodies with **any arbitrary degrees of freedom**.
* Binaries distributed via `conda <https://www.anaconda.com/download/>`_.
* Python object oriented API.
* No error message in French.
* A **cleaner code** with unit tests, more comments, no Fortran 77 and no global variables.
* OpenMP parallelization.
* **Up to 8 times faster** than Nemoh 2.0 on a single thread [#]_.
* Double precision floating numbers by default. (Can be recompiled in single precision.)
* Various input mesh formats supported via meshmagick_.
* Output in legacy Nemoh Tecplot format or NetCDF format.
* Computation of the **free surface elevation** and the **Kochin function**.
* 3D animations of the free surface elevation.

Experimental features / For advanced users
------------------------------------------

* Input via legacy ``Nemoh.cal`` files.
* Computation of the Response Amplitude Operators (RAO).
* Toeplitz matrices for faster simulations for floating bodies with local symmetries or for regular arrays of identical floating bodies.
* Hierarchical matrices with low-rank blocks for faster simulations.
* Iterative linear system solver.

Planned features (but don't expect them soon)
---------------------------------------------

* Hydrostatics.
* Removal of the irregular frequencies.


.. rubric:: Footnotes

.. [#] Version 3.6 or higher.
.. [#] However, the RAM usage of Capytaine is usually higher than the one of Nemoh 2.0. If your problem is memory bounded, Capytaine won't help (unless you can use some experimental advanced features).
