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
    
    f = np.linspace(0.1, 2.0)
    omega = 2*np.pi*f
    rho_water = 1e3
    r = 1

    sphere = Sphere(radius=r, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.center_of_mass = np.array([0, 0, 0])
    sphere.add_all_rigid_body_dofs()

    sphere.inertia_matrix = sphere.compute_rigid_body_inertia(rho=rho_water)
    sphere.hydrostatic_stiffness = sphere.compute_hydrostatic_stiffness(rho=rho_water)

    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
        'rho': rho_water,
        'water_depth': [np.infty],
        'omega': omega,
        'wave_direction': 0,
        'radiating_dof': list(sphere.dofs.keys()),
        })
    
    data = solver.fill_dataset(test_matrix, sphere_fb,
                               hydrostatics=True,
                               mesh=True,
                               wavelength=True,
                               wavenumber=True)
    
    Zi = impedance(data)



Note that we assigned the inertia and stiffness to attributes of :code:`body` called :code:`inertia_matrix` and :code:`hydrostatic_stiffness`.
These are the names expected by the :code:`fill_dataset` and :code:`impedance` functions to compute the impedance matrix.

By simple extension of incorporating the excitation transfer function response
amplitude operator (RAO)::

    from capytaine.post_pro import rao
    rao = rao(data)
