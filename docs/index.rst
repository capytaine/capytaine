Capytaine: a Python-based rewrite of Nemoh
==========================================

Capytaine_ is a boundary element method solver for the linear potential flow wave theory, written in Python and Fortran 90.

It is based on

* Nemoh_ Fortran core routines for the computation of the Green function,
* meshmagick_ for the manipulation of meshes,
* and various tools from the `Python scientific ecosystem`_.

.. _Capytaine: https://github.com/mancellin/capytaine
.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp
.. _meshmagick: https://github.com/LHEEA/meshmagick
.. _`Python scientific ecosystem`: https://scipy.org/

Features
--------

* Computation of the **added masses, radiation dampings, diffraction forces and Froude-Krylov forces** for rigid bodies or for bodies with **any arbitrary degrees of freedom**.
* Windows and Linux binaries, distributed via conda_.
* Python object oriented API.
* A **cleaner code** with unit tests, more comments, no Fortran 77 and no global variables.
* **2 to 8 times faster** than Nemoh 2.0.
* Various input mesh formats supported via meshmagick_.
* Output in legacy Nemoh Tecplot format or NetCDF format.
* Computation of the **free surface elevation** and the **Kochin function**.

.. _conda: https://www.anaconda.com/download/

Experimental features
~~~~~~~~~~~~~~~~~~~~~

* Input via legacy ``Nemoh.cal`` files.
* Faster simulations using the symmetries of the floating bodies.
* 3D animations of the free surface elevation.


Planned features
~~~~~~~~~~~~~~~~

* Fast simulation of array of identical floating bodies.
* Output in BEMIO (WEC-SIM) format.
* Hydrostatics.
* Pressure map on the body surface.

User manual
-----------

.. toctree::
   :maxdepth: 1

   user_manual/installation.rst
   user_manual/cookbook.rst

..   user_manual/tutorial.rst

Developer manual
----------------

.. toctree::
   :maxdepth: 1

   developer_manual/installation.rst
   developer_manual/code_structure.rst

Theory manual
-------------

.. toctree::
   :maxdepth: 1

   theory_manual/theory.rst
   theory_manual/bibliography.rst

License
-------

Capytaine is distributed under the terms of the GNU General Public License (GPL) v3.0. See the ``LICENSE`` file in the `code repository`_.

.. _`code repository`: https://github.com/mancellin/capytaine

This documentation is licensed under the `Creative Commons Attribution-ShareAlike 4.0 International License`_ |CCBYSA|.

.. |CCBYSA| image:: https://i.creativecommons.org/l/by-sa/4.0/80x15.png
.. _`Creative Commons Attribution-ShareAlike 4.0 International License`: http://creativecommons.org/licenses/by-sa/4.0/

.. Indices and tables
   ------------------
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
