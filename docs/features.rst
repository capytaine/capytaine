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

Below some features of the latest version (|version|) of Capytaine are listed.

Main features
-------------

* Computation of the **added masses, radiation dampings, diffraction forces and Froude-Krylov forces** for rigid bodies or for bodies with **any arbitrary degrees of freedom**.
* Windows and Linux binaries, distributed via conda_.
* Python object oriented API.
* A **cleaner code** with unit tests, more comments, no Fortran 77 and no global variables.
* **2 to 8 times faster** than Nemoh 2.0 [#]_.
* Double precision computations by default. (Single precision also possible.)
* Various input mesh formats supported via meshmagick_.
* Output in legacy Nemoh Tecplot format or NetCDF format.
* Computation of the **free surface elevation** and the **Kochin function**.

.. _conda: https://www.anaconda.com/download/

Experimental features
---------------------

* Input via legacy ``Nemoh.cal`` files.
* Faster simulations using the symmetries of the floating bodies.
* 3D animations of the free surface elevation.


Planned features
----------------

* Faster simulation of regular arrays of identical floating bodies.
* Output in BEMIO (WEC-SIM) format.
* Hydrostatics.
* Pressure map on the body surface.
* IRF computations.


.. rubric:: Footnotes

.. [#] Part of the speedup is obtained by storing intermediate computations. If your use of Nemoh is limited by the RAM usage, it might not help you.
