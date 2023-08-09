=======================
Overview of the modules
=======================


:mod:`~capytaine.matrices`
--------------------------
    Depends on: none.

    This module contains classes describing block matrices, block Toeplitz
    matrices (stored sparsely in memory) and low-rank matrices. All block
    matrices can nested to build hierarchical matrices.

    .. toctree::

       api/capytaine.matrices.block
       api/capytaine.matrices.block_toeplitz
       api/capytaine.matrices.builders
       api/capytaine.matrices.linear_solvers
       api/capytaine.matrices.low_rank

:mod:`~capytaine.meshes`
------------------------
    Depends on: none.

    This module contains classes describing the mesh of a floating body.
    Several variants, including collections of meshes and symmetric meshes are
    also available. Most of this module comes from the `meshmagick package
    <https://github.com/lheea/meshmagick>`_ by
    François Rongère.

    .. toctree::

       api/capytaine.meshes.clipper
       api/capytaine.meshes.collections
       api/capytaine.meshes.geometry
       api/capytaine.meshes.meshes
       api/capytaine.meshes.properties
       api/capytaine.meshes.quality
       api/capytaine.meshes.surface_integrals
       api/capytaine.meshes.symmetric

:mod:`~capytaine.bodies`
------------------------
    Depends on: :mod:`~capytaine.meshes`.

    This module contains a class describing a floating body, that is the reunion
    of a mesh and some degrees of freedom. It also contains some tools to
    generate new bodies of simple geometrical shapes.

    .. toctree::

       api/capytaine.bodies.bodies
       api/capytaine.bodies.predefined.spheres
       api/capytaine.bodies.predefined.cylinders
       api/capytaine.bodies.predefined.rectangles

:mod:`~capytaine.tools`
-----------------------
    Depends on: none.

    Unsorted tools.

    .. toctree::

       api/capytaine.tools.prony_decomposition

:mod:`~capytaine.green_functions`
---------------------------------
    Depends on: :mod:`~capytaine.tools`.

    This module contains the routine to evaluate the Green function.

    .. toctree::

       api/capytaine.green_functions.abstract_green_function
       api/capytaine.green_functions.delhommeau

:mod:`~capytaine.bem`
----------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.matrices`, :mod:`~capytaine.green_functions`, :mod:`io.xarray <capytaine.io.xarray>`.

    The module is the core of the code. It contains the routines to assemble the
    matrices and solve the BEM problem.

    .. toctree::

       api/capytaine.bem.airy_waves
       api/capytaine.bem.engines
       api/capytaine.bem.solver
       api/capytaine.bem.problems_and_results

:mod:`io.xarray <capytaine.io.xarray>`
---------------------------------------
    Depends on: :mod:`~capytaine.bem`.

    This submodule contains the code used to read and write the :code:`xarray`
    datasets that are the standard output of the code. It is interlaced with
    :code:`capytaine.bem` and might be merged with it in the future.

    .. toctree::

       api/capytaine.io.xarray

:mod:`~capytaine.post_pro`
--------------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`, :mod:`io.xarray <capytaine.io.xarray>`.

    This module contains various tools for the post-processing of the results of
    the BEM problem.

    .. toctree::

       api/capytaine.post_pro.free_surfaces
       api/capytaine.post_pro.kochin
       api/capytaine.post_pro.impedance
       api/capytaine.post_pro.rao

:mod:`ui.vtk <capytaine.ui.vtk>`
--------------------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies` and :mod:`~capytaine.post_pro`.

    This module contains the calls to VTK used for the 3D display of the meshes
    and the animation of the free surface.

    .. toctree::

       api/capytaine.ui.vtk.animation
       api/capytaine.ui.vtk.body_viewer
       api/capytaine.ui.vtk.helpers
       api/capytaine.ui.vtk.mesh_viewer

:mod:`~capytaine.io`
--------------------
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`.

    This module contains various tools for inputs and outputs of the code.

    .. toctree::

       api/capytaine.io.legacy
       api/capytaine.io.mesh_loaders
       api/capytaine.io.mesh_writers

:mod:`ui.cli <capytaine.ui.cli>`
---------------------------------

Depends on: :mod:`~capytaine.io`, :mod:`~capytaine.bem`.

    This module contains the command-line interface of the code.

    .. toctree::

       api/capytaine.ui.cli
