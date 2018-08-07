==============
Code structure
==============

Core of the code
----------------

The core of the code is ``Nemoh.py`` and the Fortran files in ``NemohCore``. 

The Fortran files are mostly used to evaluate the Green function and integrate it to build the influence matrices. This is split into two files: ``Green_Rankine.f90`` deals with the :math:`1/r` part of the Green function whereas ``Green_wave.f90`` deals with the wave part of the Green function.

In ``Nemoh.py`` the linear system is solved and the added mass, radiation damping or diffraction are deduced.
In this file, the decorators :code:`@use_symmetries` and :code:`@lru_cache` are optimizations of the code.
For debugging, they can be removed and the code should run the same (but more slowly).

The files ``bodies.py``, ``problems.py`` and ``results.py`` are supporting modules.
The first one contains the ``FloatingBody`` class which is basically the union of a ``Mesh`` from the ``meshmagick`` library and some degrees of freedom.
The latter two contain classes storing the linear potential flow problems and their solutions.

Some simple floating bodies with geometrical shapes can be generated using the submodule ``geometric_bodies``.

More advanced features
----------------------

The files ``symmetries.py`` and ``Toeplitz_matrices.py`` are used to speed up the computations for symmetric bodies (see [AD18]_).
The submodule ``tools.vtk`` contains some routines to display the mesh and the free surface in 3D using VTK.

User interface
--------------

The file ``cli.py`` contains an experimental command-line interface.
The module ``tools.import_export`` contains also some helper functions to interface with Nemoh 2 format.

