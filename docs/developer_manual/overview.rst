=======================
Overview of the modules
=======================


:mod:`~capytaine.matrices`
    
:mod:`~capytaine.meshes`

:mod:`~capytaine.bodies` 
    Depends on: :mod:`~capytaine.meshes`.

:mod:`~capytaine.bem` 
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.matrices`, :mod:`io.xarray <capytaine.io.xarray>`.

:mod:`io.xarray <capytaine.io.xarray>` 
    Depends on: :mod:`~capytaine.bem`.

:mod:`~capytaine.post_pro` 
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`, :mod:`io.xarray <capytaine.io.xarray>`.

:mod:`ui.vtk <capytaine.ui.vtk>` 
    Depends on: :mod:`~capytaine.meshes`.

:mod:`~capytaine.io`
    Depends on: :mod:`~capytaine.meshes`, :mod:`~capytaine.bodies`, :mod:`~capytaine.bem`.

:mod:`ui.cli <capytaine.ui.cli>` 
    Depends on: :mod:`~capytaine.io`, :mod:`~capytaine.bem`.
