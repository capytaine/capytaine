===============
Post-processing
===============

Free surface elevation
----------------------

To compute the free surface elevation, let us first initialize a
:class:`~capytaine.post_pro.free_surfaces.FreeSurface` object::

    from capytaine import FreeSurface
    fs = FreeSurface(x_range=(-10, 10), nx=10, y_range=(-5, 5), ny=10)

The above code generates a regular free surface mesh of :math:`10 \times 10`
cells. This object can be used by the solver with a results object to get the
free surface elevation::

    fs_elevation = solver.get_free_surface_elevation(result, free_surface)

The output is a numpy array storing the free surface elevation in frequency
domain as a complex number at each point of the free surface (in the present
example an array of shape :code:`(10, 10)`).

The result object should have been computed with the option
:code:`keep_details=True`. The solver does not need to be the one that computed
the result object.

The undisturbed incoming waves (Airy waves) can be computed as follow::

    incoming_waves = fs.incoming_waves(DiffractionProblem(omega=1.0, angle=pi/2))

See the examples in the :doc:`cookbook` for usage in a 3D animation.

Impedance and RAO
-----------------

The intrinsic impedance can be computed based on the hydrodynamics, 
hydrostatics, and inertial properties::

    import numpy as np
    import xarray as xr
    from capytaine import BEMSolver
    from capytaine.bodies.predefined.spheres import Sphere
    from capytaine.post_pro import impedance
    
    sphere = Sphere(radius=r, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_all_rigid_body_dofs()
    
    f = np.linspace(0.1, 2.0)
    omega = 2*np.pi*f
    rho_water = 1e3
    r = 1
    m = 1.866e+03
    
    M = np.array([
           [ m,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00],
           [ 0.000e+00,  1.866e+03,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00],
           [ 0.000e+00,  0.000e+00,  1.866e+03,  0.000e+00,  0.000e+00,  0.000e+00],
           [ 0.000e+00,  0.000e+00,  0.000e+00,  4.469e+02,  9.676e-31, -2.757e-14],
           [ 0.000e+00,  0.000e+00,  0.000e+00,  9.676e-31,  4.469e+02,  3.645e-15],
           [ 0.000e+00,  0.000e+00,  0.000e+00, -2.757e-14,  3.645e-15,  6.816e+02]])
    
    kHS = np.array([
        [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
        [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ],
        [    0.   ,     0.   , 29430.   ,     0.   ,     0.   ,     0.   ],
        [    0.   ,     0.   ,     0.   ,   328.573,     0.   ,     0.   ],
        [    0.   ,     0.   ,     0.   ,     0.   ,   328.573,     0.   ],
        [    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ]])
    
    sphere.mass = sphere.add_dofs_labels_to_matrix(M)
    sphere.hydrostatic_stiffness = sphere.add_dofs_labels_to_matrix(kHS)
    solver = BEMSolver()
    
    test_matrix = xr.Dataset(coords={
        'rho': rho_water,                         
        'water_depth': [np.infty],          
        'omega': omega,
        'wave_direction': 0,
        'radiating_dof': list(sphere_fb.dofs.keys()),
        })
    
    data = solver.fill_dataset(test_matrix, [sphere_fb],
                               hydrostatics=True,
                               mesh=True,
                               wavelength=True,
                               wavenumber=True)
    
    Zi = impedance(data)

By simple extension of incorporating the excitation transfer function response 
amplitude operator (RAO)::

    from capytaine.post_pro import rao
    rao = rao(data)
