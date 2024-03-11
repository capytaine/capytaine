===============
Post-processing
===============

Pressure, velocity, free surface elevation
------------------------------------------

Once the problem has been solved, several fields of interest can be computed at post-processing:

+-----------------------------------------------------------+------------------------------------------------------+
| Code                                                      | Description                                          |
+===========================================================+======================================================+
| ``solver.compute_potential(points, result)``              | The velocity potential :math:`\phi(x, y, z)`         |
+-----------------------------------------------------------+------------------------------------------------------+
| ``solver.compute_pressure(points, result)``               | The pressure in the fluid :math:`p(x, y, z)`         |
+-----------------------------------------------------------+------------------------------------------------------+
| ``solver.compute_velocity(points, result)``               | The velocity of the fluid :math:`u(x, y, z)`         |
+-----------------------------------------------------------+------------------------------------------------------+
| ``solver.compute_free_surface_elevation(points, result)`` | The elevation of the free surface :math:`\eta(x, y)` |
+-----------------------------------------------------------+------------------------------------------------------+

All the methods listed above work in the same way: they require the :class:`~capytaine.bem.problems_and_results.LinearPotentialFlowResult` object containing the required data about the solved problem and some points at which the field should be evaluated.

The result object should have been computed with the indirect method (on by default) and the option :code:`keep_details=True` (on by default for :meth:`~capytaine.bem.solver.BEMSolver.solve`, off by default for :meth:`~capytaine.bem.solver.BEMSolver.solve_all`).
The solver does not need to be the one that computed the result object.

.. note::
    The functions in the :mod:`~capytaine.bem.airy_waves`, used to compute the same magnitudes for an undisturbed incoming wave field, have the same structure.

The point(s) can be given in several ways:

- Either a single point, given as a list, a tuple, or an 1d-array::

    solver.compute_potential([3.0, -2.0, -5.0], results)

- or a list of points, given as a list of lists, or a list of tuples, or a 2d-array::

    solver.compute_potential([[3.0, -2.0, -5.0], [4.0, 5.0, -2.0]], results)

- or the return of a call to ``meshgrid``::

    points = np.meshgrid(np.linspace(-2.0, 2.0, 10), np.linspace(-3.0, 3.0, 20), np.linspace(-2.0, 0.0, 30))
    solver.compute_potential(points, results)

- or a mesh, in which case the centers of the faces of the mesh are used::

    solver.compute_potential(mesh, results)

- or a floating body, in which case the corresponding mesh will be used::

    solver.compute_potential(body, results)

- or a :class:`~capytaine.post_pro.free_surfaces.FreeSurface` object, although the use of this object is not recommended unless you are preparing a 3D animation with the Capytaine's VTK viewer which still require this object at the moment::

    fs = cpt.FreeSurface(x_range=(-10, 10), y_range=(-10, 10))
    solver.compute_potential(fs, results)

The returned values is an array of shape matching the shape of the input points.

.. warning::
   There is a single case in which passing a mesh is not equivalent to a list of point: if you want the compute the velocity on the hull of the floating body. In this case, you should give the same mesh object that has been used for the resolution::

        solver.compute_velocity(result.body.mesh, result)

   Other Python objects might return incorrect values or errors.

For potential, pressure and velocity, 3 coordinates :math:`(x, y, z)` are expected for each points.
For the free surface elevation, 2 coordinates :math:`(x, y)` are sufficient.

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
        'water_depth': [np.inf],
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
