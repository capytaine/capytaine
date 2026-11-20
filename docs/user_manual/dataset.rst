=====================
Collection in dataset
=====================

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
contains the potential and pressure fields on each face of the mesh used for
the resolution (that is, including the lid faces, if a lid has been defined
on the body to remove irregular frequencies)::

  print(result.potential)
  # [-1.72534485e-03+0.14128629j -3.98932611e-03+0.33387497j
  #  -6.54135830e-03+0.54208628j ... -2.09853806e+00-0.56106653j
  #  -2.18441640e+00-0.58402701j -2.22777755e+00-0.59562008j]

  print(result.pressure)
  # [-1.72534485e-03+0.14128629j -3.98932611e-03+0.33387497j
  #  -6.54135830e-03+0.54208628j ... -2.09853806e+00-0.56106653j
  #  -2.18441640e+00-0.58402701j -2.22777755e+00-0.59562008j]

These magnitudes are stored in an one-dimensional array as long as the number
of faces of the mesh, and stored in the same order as the faces in the
``mesh_including_lid`` object. In other words, ``result.pressure[3]``
contains the pressure on the face of center
``body.mesh_including_lid.faces_centers[3]``.

Since the pressure on the lid is not physically meaningful, the result object
also provides ``result.pressure_on_hull``, which is ``result.pressure``
restricted to the faces of ``body.mesh`` (that is, excluding the lid). When no
lid has been defined, ``result.pressure_on_hull`` and ``result.pressure`` are
the same.

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

- :code:`hydrostatics` (default: :code:`True`): tries to compute the
  hydrostatic properties of the ``FloatingBody`` and add them to the dataset.

- :code:`mesh` (default: :code:`False`): stores in the dataset the full mesh for
  the hull and the lid, as well as the dofs.

  - The hull mesh is stored using the ``hull_face`` dimension in the arrays :code:`mesh_vertices` (shape ``(hull_face, 4, 3)``) and :code:`mesh_faces_center` (shape: ``(hull_face, 3)``).

  - If the body has a lid, it is stored similarly using the ``lid_face`` dimension in the arrays :code:`lid_mesh_vertices` (shape ``(lid_face, 4, 3)``) and :code:`lid_mesh_faces_center` (shape: ``(lid_face, 3)``).

  - If the body has dofs, the array ``dof_motions`` (shape: ``(influenced_dof, hull_face, 3)``) stores the displacement at each face center due to each dof.

  - If the dofs have a Jacobian defined (as do rigid body dofs or :class:`~capytaine.bodies.dofs.CustomDof` when setting the ``gradient_of_motion`` parameter), it is stored also in an array of shape ``(influenced_dof, hull_face, 3, 3)``.

  - The name of the quadrature scheme used for the mesh is stored as a string in the :code:`quadrature_method` variable.

.. note:: If the body's mesh takes advantage of a planar or rotational
          symmetry (e.g. a :class:`~capytaine.meshes.symmetric_meshes.ReflectionSymmetricMesh`),
          the symmetry itself is **not** preserved by :code:`mesh=True`, as of
          the current version of Capytaine: :code:`mesh_vertices`,
          :code:`mesh_faces_center` and :code:`dof_motions` always contain the
          full, explicit list of faces (as if :meth:`~capytaine.meshes.AbstractMesh.merged`
          had been called first), with nothing in the dataset indicating that
          the original mesh was symmetric. This only affects the exported
          mesh/dofs data itself; the hydrodynamic coefficients are unaffected
          and are still computed by taking advantage of the symmetry
          internally.

- In addition, :code:`fill_dataset` accepts a :code:`keep_details` keyword
  argument (default: :code:`False`), forwarded to :meth:`~capytaine.bem.solver.BEMSolver.solve_all`.
  When set to :code:`True`, the potential, pressure and sources computed for
  each problem is stored in the dataset:

   - the variables :code:`pressure_diffraction`, :code:`pressure_Froude_Krylov` and :code:`pressure_radiation` contain the pressure on the hull (with a :code:`hull_face` dimension, in the same order as the faces of :code:`body.mesh`).
     Pressure on the lid has no physical meaning and is not stored in the dataset.

   - if the problems were solved with the indirect solver, the dataset includes the variables ``sources_diffraction`` and ``sources_radiation`` (on the hull, with dimension ``hull_face``) and ``sources_lid_diffraction`` and ``sources_lid_radiation`` (on the lid, if any, with dimension ``lid_face``)

   - with the direct method, the potential is stored instead, as ``potential_diffraction``, ``potential_radiation``, ``potential_lid_diffraction`` and ``potential_lid_radiation``. Results of the direct and indirect methods solved with ``keep_details=True`` cannot be assembled in the same dataset (they can with ``keep_details=False``).

Keeping these details increases the memory usage of the dataset, so this
option is disabled by default.

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

The results can also be collected by :func:`~capytaine.io.xarray.assemble_matrices`, which returns some matrices of :func:`~capytaine.io.xarray.assemble_dataset` as numpy arrays stripped of their metadata.
This function is meant to be used for teaching, to assemble the matrices without getting the students in contact with ``xarray``.
