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
    The functions in the :mod:`~capytaine.bem.airy_waves` module, used to compute the same magnitudes for an undisturbed incoming wave field, have the same structure.

The point(s) can be given in several ways:

- Either a single point, given as a list, a tuple, or an 1d-array::

    solver.compute_potential([3.0, -2.0, -5.0], result)

- or a list of points, given as a list of lists, or a list of tuples, or a 2d-array::

    solver.compute_potential([[3.0, -2.0, -5.0], [4.0, 5.0, -2.0]], result)

- or the return of a call to ``meshgrid``::

    points = np.meshgrid(np.linspace(-2.0, 2.0, 10), np.linspace(-3.0, 3.0, 20), np.linspace(-2.0, 0.0, 30))
    solver.compute_potential(points, result)

- or a mesh, in which case the centers of the faces of the mesh are used::

    solver.compute_potential(mesh, result)

in particular, the following can be used to compute the value on the mesh used for the resolution::

    solver.compute_potential(result.body.mesh, result)

- or a floating body, in which case the corresponding mesh will be used::

    solver.compute_potential(body, result)

The returned values is an array of shape matching the shape of the input points.

For potential, pressure and velocity, 3 coordinates :math:`(x, y, z)` are expected for each points.
For the free surface elevation, 2 coordinates :math:`(x, y)` are sufficient.

.. note::
   In the limit cases of zero and infinite frequencies, these magnitudes can also be computed.
   Strictly speaking, their value is 0 (resp. infinity) at these frequencies.
   Nonetheless, the following magnitudes :math:`\phi/\omega`, :math:`u/\omega`, :math:`\eta/\omega^2` and :math:`p/\omega^2` are well defined and have a finite value.
   They can be accessed in Capytaine as in the following example::

       [...]
       pb = cpt.RadiationProblem(body=body, omega=0, radiating_dof="Heave")
       res = solver.solve(pb, keep_details=True)
       pressure_over_omega_2 = solver.compute_pressure(points, res)/(pb.omega*pb.omega)

.. warning::
   The velocity cannot be reliably reconstructed at points lying on a face of the mesh that is on the free surface
   (typically a face of a lid used for irregular frequencies removal), except at the center of the face (e.g. when
   passing the lid mesh itself as points).
   The singular part of the Green function on the free surface is only integrated specifically when the point is at the
   center of the face. At other points of the face, the velocity is inaccurate at finite frequency and wrong at
   infinite frequency.
   Even at the center of the face, the face is approximated by a disk of the same area for this integral.
   Such points are inside the body, that is outside of the fluid domain, so this is not expected to be a limitation in
   practice. Points on the free surface outside of the body are not affected.


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

Kochin functions
----------------

The function :func:`~capytaine.post_pro.kochin.compute_kochin` can compute the
Kochin function of the waves for a
:class:`~capytaine.bem.problems_and_results.LinearPotentialFlowResult`
containing the source distribution (that is, if it has been solved with the
indirect method)::

    from capytaine.post_pro.kochin import compute_kochin
    thetas = np.linspace(0, 2*np.pi, 20)
    kochin_values = compute_kochin(result, thetas)

Alternatively, if the ``test_matrix`` of
:func:`~capytaine.bem.solver.BEMSolver.fill_dataset` contains a coordinate called
``theta``, it is used to compute the Kochin function of all solved problems
(see cookbook example).

Mean drift forces
-----------------

Near-field
~~~~~~~~~~

The mean drift forces can be computed with :func:`~capytaine.post_pro.mean_drift_force.near_field_mean_drift_force`.
Here is an example::

    import numpy as np
    import xarray as xr
    import capytaine as cpt

    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    solver = cpt.BEMSolver(method='indirect')
    test_matrix = xr.Dataset(coords={
                'wavenumber': [1.16, 1.24, 1.39],
                'wave_direction': [0.0],
                'radiating_dof': list(body.dofs.keys())
            })
    dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True, keep_details=True)
    first_order_motion = cpt.post_pro.rao(dataset)

    from capytaine.post_pro.mean_drift_force import near_field_mean_drift_force
    nf_mdf = near_field_mean_drift_force(first_order_motion, dataset, body, solver=solver, output_pressure=False)

The ``dataset`` should have been computed with the ``'indirect'`` method and ``keep_details=True``, such that the source distribution is stored in the dataset allowing first-order quantities to be used in the computation of the mean drift force.
The force is computed for the frequencies and wave directions common to ``first_order_motion_amplitude`` and ``dataset``, and for the dofs common to them and ``body`` (a warning is logged when some are left out).
This function returns a dataset with a ``near_field_mean_drift_force`` variable of shape ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_dofs)``.

The implementation supports single rigid body and multiple bodies.
The expression used for the computation of the mean drift force are generic enough to be used with generalized degrees of freedom, but this use case has not been validated.

Optional parameter ``output_pressure`` allows to return the mean drift pressure field on the hull as ``pressure_hull_mean_drift``, of shape ``(nb_freq, nb_wave_direction, nb_wave_direction, nb_hull_faces)``, as well as the first-order fields on the waterline from which the waterline term of the force is computed: ``pressure_waterline_total`` and ``relative_elevation_waterline``, of shape ``(nb_freq, nb_wave_direction, nb_waterline_segments)``, and the waterline segments ``mesh_waterline``.

Far-field
~~~~~~~~~

The horizontal mean drift forces can be computed from the Kochin function using
:func:`~capytaine.post_pro.mean_drift_force.far_field_mean_drift_force`,
here is an example::

    import numpy as np
    import xarray as xr
    import capytaine as cpt

    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))
    solver = cpt.BEMSolver(method='indirect')
    test_matrix = xr.Dataset(coords={
                'wavenumber': [1.16, 1.24, 1.39],
                'wave_direction': [0.0],
                'theta': np.linspace(-np.pi/8, 2*np.pi, 100),
                'radiating_dof': list(body.dofs.keys())
            })
    dataset = solver.fill_dataset(test_matrix, body, hydrostatics=True)
    first_order_motion = cpt.post_pro.rao(dataset)

    from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force
    ff_mdf = far_field_mean_drift_force(first_order_motion, dataset)

Here we added ``theta`` to the ``test_matrix`` such that the Kochin function is
automatically included in the dataset.
The ``'indirect'`` method (default) is required for this.
The range of ``theta`` should cover both :math:`[0, 2\pi]` and some
neighborhood of the wave direction of interest (here :math:`[-\pi/8, \pi/8]`).

The output is a dataset with variables ``drift_force_surge``, ``drift_force_sway`` and ``drift_force_yaw``.
They can be stacked into a single variable with dimension ``influenced_dof``
like the other forces in Capytaine with
:func:`~capytaine.post_pro.mean_drift_force._merge_far_field_mean_drift_variables`.
