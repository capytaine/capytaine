=======================
Overview of the modules
=======================

:mod:`~capytaine.tools`
-----------------------
    Depends on: none.

    Miscellaneous independent tools, used by most of the other modules.
    Especially :mod:`~capytaine.tools.optional_imports` and
    :mod:`~capytaine.tools.deprecation_handling` are used all over the place
    without being explicitly referenced as a dependency in the present page.

    .. toctree::
       :maxdepth: 1

       api/capytaine.tools


:mod:`~capytaine.matrices`
--------------------------
    Depends on: none.

    This module contains data structures to represent matrices as block
    matrices, block Toeplitz matrices (stored sparsely in memory) and low-rank
    matrices. All block matrices can nested to build hierarchical matrices.

    These matrices representations are advanced features that are not used by
    default by Capytaine. In the future, this module could be extracted as an
    indendant library to make Capytaine a bit leaner.

    .. toctree::
       :maxdepth: 1

       api/capytaine.matrices


:mod:`~capytaine.meshes`
------------------------
    Depends on: none.

    This module contains classes describing the mesh of a floating body.
    Several variants, including collections of meshes and symmetric meshes are
    also available. Most of this module comes from the `meshmagick package
    <https://github.com/lheea/meshmagick>`_ by François Rongère.  It also
    contains some tools to generate meshes of simple geometrical shapes.


    .. toctree::
       :maxdepth: 1

       api/capytaine.meshes


:mod:`~capytaine.bodies`
------------------------
    Depends on: :mod:`~capytaine.meshes`.

    This module contains a class describing a floating body, that is the reunion
    of a mesh and some degrees of freedom, as well as some more optional data.

    .. toctree::
       :maxdepth: 1

       api/capytaine.bodies


:mod:`~capytaine.green_functions`
---------------------------------
    Depends on: :mod:`~capytaine.tools`.

    This module contains the routine to evaluate the Green function.

    .. toctree::
       :maxdepth: 1

       api/capytaine.green_functions


:mod:`~capytaine.bem`
----------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.matrices`, :mod:`~capytaine.green_functions`, :mod:`io.xarray <capytaine.io.xarray>` and :mod:`~capytaine.tools`.

    The module is the core of the code. It contains the routines to assemble the
    matrices and solve the BEM problem.

    .. toctree::
       :maxdepth: 1

       api/capytaine.bem


:mod:`io.xarray <capytaine.io.xarray>`
---------------------------------------
    Depends on: :mod:`~capytaine.bem`.

    This submodule contains the code used to read and write the :code:`xarray`
    datasets that are the standard output of the code. It is interlaced with
    :code:`capytaine.bem` and might be merged with it in the future.

    .. toctree::
       :maxdepth: 1

       api/capytaine.io.xarray


:mod:`~capytaine.post_pro`
--------------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`, :mod:`io.xarray <capytaine.io.xarray>`.

    This module contains various tools for the post-processing of the results of
    the BEM problem.

    .. toctree::
       :maxdepth: 1

       api/capytaine.post_pro


:mod:`~capytaine.io`
--------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`.

    This module contains various tools for inputs and outputs of the code.

    .. toctree::
       :maxdepth: 1

       api/capytaine.io


:mod:`ui <capytaine.ui>`
--------------------------------
    Depends on most of the other modules.

    This modules contains the code handling the user interfaces:
    the display of the outputs in the terminal using Rich, the command line interface
    and the 3D visualisations of the meshes with VTK.

    .. toctree::
       :maxdepth: 1

       api/capytaine.ui
