=========================================
Post-processing and collection in dataset
=========================================

Outputs stored in LinearPotentialFlowResult
-------------------------------------------

Each ``LinearPotentialFlowResult`` object contains all the data defined for the
corresponding problem::

  problem = cpt.DiffractionProblem(body=body, period=1.0, wave_direction=3.14)
  solver = cpt.BEMSolver()
  result = solver.solve(problem, keep_details=True)

  print(result.period)
  # 1.0

They also contain a supplementary attribute named ``forces``, containing the
computed forces on the body::

  print(result.forces)
  # {'Surge': (23941.728129534735+9487.048661471363j), 'Sway': (-2.010438038269058e-09-2.8862814360763878e-09j), ...}

The force is stored as Python dictionary associating a complex value in SI
units to the name of the influenced degree of freedom.
In other words, in the example code above, the :math:`x`-component of the force
vector has magnitude 23941 Newton at :math:`t =0`, while the
:math:`y`-component of the force vector is negligible at all time.

For radiation problems, the result object also contain ``added_mass`` and
``radiation_damping`` attributes, with the same shape of a Python dictionary.

It the solver was called with ``keep_details=True``, the result object also
contains the potential and pressure fields on each face of the mesh of the hull::

  print(result.potential)
  # [-1.72534485e-03+0.14128629j -3.98932611e-03+0.33387497j
  #  -6.54135830e-03+0.54208628j ... -2.09853806e+00-0.56106653j
  #  -2.18441640e+00-0.58402701j -2.22777755e+00-0.59562008j]

  print(result.pressure)
  # [-1.72534485e-03+0.14128629j -3.98932611e-03+0.33387497j
  #  -6.54135830e-03+0.54208628j ... -2.09853806e+00-0.56106653j
  #  -2.18441640e+00-0.58402701j -2.22777755e+00-0.59562008j]

These magnitudes are stored in an one-dimensional array as long as the number
of faces of the mesh, and stored in the same order as the faces in the ``Mesh``
object. In other words, ``result.pressure[3]`` contains the pressure on the
face of center ``body.mesh.faces_centers[3]``.

Recall that the potential and the pressure are related by :math:`p = j \rho
\omega \Phi`, where :math:`\rho` is the fluid density and :math:`\omega` is the
angular frequency.

The potential, the pressure and the velocity in the rest of the fluid can be
computed in post-processing as described below.

When using the indirect boundary integral equation (by default), the result
object with details also contains the distribution of sources :math:`\sigma` on
the hull, under the same form as the potential above::

  print(result.sources)
  # [  0.0281344 -0.35719578j   0.06715056-1.34127138j
  #    0.10891061-2.26921669j ... -13.58069557+2.13219904j
  #  -14.13645749+2.21945488j -14.41706935+2.26351156j]


Building a dataframe or a dataset from LinearPotentialFlowResult
----------------------------------------------------------------

If you have a list of :code:`LinearPotentialFlowResult`, you can assemble
them in a `pandas <https://pandas.pydata.org/>`_ DataFrame or a `xarray
<https://docs.xarray.dev>`_ dataset for a more convenient post-processing. Use
the following syntax::

   dataframe = cpt.assemble_dataframe(list_of_results)

or::

   dataset = cpt.assemble_dataset(list_of_results)

If you gave a test matrix to the :code:`BEMSolver.fill_dataset` method, the
output will directly be an xarray dataset.

Errored resolutions stored as a
:class:`~capytaine.bem.problems_and_results.FailedDiffractionResult` or a
:class:`~capytaine.bem.problems_and_results.FailedRadiationResult` will appear
as `NaN` in the dataset.

Both :code:`assemble_dataset` and :code:`fill_dataset` accept some optional keyword
arguments to store more information in the dataset:

- :code:`wavenumber`, :code:`wavelength`, :code:`period`, :code:`omega`,
  (default: all `True`): control whether which of the representations of the
  wave frequency are stored in the dataset. At least one should be included, by
  default they all are.
- :code:`mesh` (default: :code:`False`): add some information about the mesh in
  the dataset (number of faces, quadrature method).
- :code:`hydrostatics` (default: :code:`True`): if hydrostatics data are
  available in the :code:`FloatingBody`, they are added to the dataset.

.. note:: The code does its best to keep the degrees of freedom in the same
          order as they have been provided by the user, but there is no
          guarantee that this will always be the case.
          Nonetheless, in all cases the labels will always match the data.
          So selecting a value by the name of the dof will always return the right one::

              data.sel(radiating_dof="Heave", influenced_dof="Heave")

          You can also manually reorder the dofs with the following syntax::

              sorted_dofs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]
              print(data.sel(radiating_dof=sorted_dofs, influenced_dof=sorted_dofs))

.. note:: Datasets created with :code:`assemble_dataset` only include data on
          cases with a free surface.
          Cases without a free surface (:code:`free_surface=inf`) are ignored.

The results can also be collected by :func:`~capytaine.io.xarray.assemble_matrices`, which returns the matrices of :func:`~capytaine.io.xarray.assemble_dataset` as numpy arrays stripped of their metadata.
This function is meant to be used for teaching, to assemble the matrices without getting the students in contact with ``xarray``.


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

- or a :class:`~capytaine.post_pro.free_surfaces.FreeSurface` object, although the use of this object is not recommended unless you are preparing a 3D animation with the Capytaine's VTK viewer which still require this object at the moment::

    fs = cpt.FreeSurface(x_range=(-10, 10), y_range=(-10, 10))
    solver.compute_potential(fs, result)

The returned values is an array of shape matching the shape of the input points.

.. warning::
   There is a single case in which passing a mesh is not equivalent to a list of point: if you want the compute the velocity on the hull of the floating body. In this case, you should give the same mesh object that has been used for the resolution::

        solver.compute_velocity(result.body.mesh, result)

   Other Python objects might return incorrect values or errors.

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
