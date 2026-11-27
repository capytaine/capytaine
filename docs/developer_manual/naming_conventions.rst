==================
Naming conventions
==================

This page lists some conventions followed by the public interface of Capytaine.
New functions should follow them. Some older functions do not, but are kept for backward compatibility.

Module names
------------

Module names are lowercase nouns.
Subpackages are named after the domain they cover,
e.g. ``meshes``, ``bodies``, ``green_functions``, ``bem``, ``post_pro``.

A module defining a family of classes is named after them in plural form,
e.g. ``meshes/symmetric_meshes.py``, ``bem/linear_solvers.py``,
with their abstract base classes in ``abstract_<classes>.py``, e.g. ``bodies/abstract_bodies.py``.

Function names
--------------

Functions and methods performing an action start with a verb describing what they do:

- ``compute_<quantity>`` for the computation of a physical quantity,
  e.g. :func:`~capytaine.post_pro.kochin.compute_kochin`,
  :func:`~capytaine.post_pro.mean_drift_force.compute_near_field_mean_drift_force`,
  :meth:`~capytaine.bem.solver.BEMSolver.compute_pressure`,
  :meth:`~capytaine.bodies.bodies.FloatingBody.compute_hydrostatics`.
- ``solve`` and ``fill_dataset`` for running the BEM solver.
- ``assemble_<something>`` for gathering results into a dataset or a data frame,
  e.g. :func:`~capytaine.io.xarray.assemble_dataset`.
  (Exception: the methods ``assemble_regular_array`` and ``assemble_arbitrary_array`` build an array of bodies.)
- ``export_<something>``, ``export_to_<format>``, ``load_<something>``, ``save_<something>`` for input and output.
- ``mesh_<shape>`` for the constructors of predefined meshes, e.g. :func:`~capytaine.meshes.predefined.spheres.mesh_sphere`.
- ``build_<something>``, ``evaluate_<something>`` and ``integrate_<something>`` for lower-level methods
  of matrices, dofs and bodies.

A noun as function name (such as ``rao``) takes a name that the user would like to use for the result::

    rao = rao(dataset)  # The variable may shadow the function

Some post-processing functions of the ``post_pro`` module, such as :func:`~capytaine.post_pro.rao.rao` and
:func:`~capytaine.post_pro.impedance.impedance`, do not follow this convention for historical reasons.


Methods and properties of objects
---------------------------------

Meshes, bodies and other objects are not modified in place: their methods return a new object.

(Exception: ``add_translation_dof``, ``add_rotation_dof``, ``add_all_rigid_body_dofs``, which are kept for internal
use and backward compatibility, but the recommended user interface is :func:`~capytaine.bodies.dofs.rigid_body_dofs`.)

Besides the verbs above, the following patterns are used:

- A past participle for a method returning a transformed copy of the object,
  e.g. ``translated``, ``rotated_x``, ``mirrored``, ``clipped``, ``merged``::

      new_mesh = mesh.translated((1.0, 0.0, 0.0)).mirrored("xOz")

- ``with_<something>`` and ``without_<something>`` for a copy with a modified attribute,
  e.g. ``with_only_dofs``, ``with_quadrature``, ``without_metadata``.
- ``from_<something>`` for the class methods building an object from another representation,
  e.g. ``Mesh.from_list_of_faces``.

A method returning a new object of another nature, derived from the object, follows the verb convention of the
functions, e.g. ``generate_lid`` returns a lid mesh for a hull mesh.

Some methods returning a modified copy do not follow these patterns for historical reasons, e.g. ``immersed_part``.

Properties are nouns, e.g. ``faces_centers`` or ``mesh_including_lid``.
Properties of the elements of a mesh are named ``<elements>_<quantity>``,
e.g. ``faces_centers``, ``faces_normals``, ``waterline_edges_centers``, ``waterline_edges_normals``.

The number of something is ``nb_<something>``, e.g. ``nb_faces``, ``nb_dofs``
(except ``n_jobs`` and ``n_threads``, following the convention of the libraries they are passed to).


Order of the arguments
----------------------

1. The main data comes first: the dataset, or otherwise the result, the problem or the body.
2. Then the other required objects it is combined with, such as the body.
3. Then the options and the secondary data, as keyword-only arguments (after ``*``).

For instance::

    solver.fill_dataset(dataset, body, *, method=None, ...)
    compute_far_field_mean_drift_force(dataset, *, rao)
    compute_near_field_mean_drift_force(dataset, body, *, rao, solver=None, output_pressure=False)

Keyword-only arguments allow adding new arguments, or giving a default value to an argument, without changing the
existing calls. They also keep the signatures of related functions consistent with each other, such as the two mean
drift force functions above.

Exception: functions evaluating a field at some points take the points first:
:meth:`~capytaine.bem.solver.BEMSolver.compute_potential` ``(points, result)``
and :func:`~capytaine.bem.airy_waves.airy_waves_pressure` ``(points, problem)``.


Dimensions and variables of the datasets
----------------------------------------

The dimensions of xarray objects are named after what indexes them, in singular form:

- The frequency is one of ``omega``, ``freq``, ``period``, ``wavelength`` or ``wavenumber``,
  as given by the user in the test matrix or in the problem objects.
- ``wave_direction``, and ``other_wave_direction`` for the second wave of the second-order quantities depending on a
  pair of waves (such as the mean drift force). This way, the first-order quantities and the second-order quantities share
  the coordinate ``wave_direction``, and a first-order quantity is moved to the second wave with
  ``.rename({"wave_direction": "other_wave_direction"})``.
- ``radiating_dof`` for the dof of the motion and ``influenced_dof`` for the dof of the force.
- ``hull_face`` and ``lid_face`` for the fields on the faces of the meshes,
  ``waterline_edge`` for the fields on the edges of the waterline.
- ``vertices_of_face`` and ``vertices_of_edge`` for the vertices of a face or of an edge of the meshes.
- ``space_coordinate`` (``x``, ``y``, ``z``) for the coordinates of a point,
  ``vector_component`` (``x``, ``y``, ``z``) for the components of a vector (e.g. a motion, a normal, a velocity),
  and ``derivative_direction`` (``x``, ``y``, ``z``) for the direction of the derivative in a Jacobian matrix.
- ``body`` for the bodies of a multibody (e.g. their centers of buoyancy).
- ``theta`` for the direction of the Kochin function.
- ``complex`` (``re``, ``im``) for the real and imaginary parts, when saving a dataset with
  :func:`~capytaine.io.xarray.separate_complex_values`.

Variables of the datasets are usually named after the pattern ``quantity_[location]_[description]``,
e.g. ``pressure_hull_diffraction``, ``sources_lid_radiation``, ``kochin_radiation``, ``relative_elevation_waterline``.
The geometry of the meshes is stored as ``mesh_<location>_<property of the mesh>``,
e.g. ``mesh_hull_vertices``, ``mesh_lid_faces_centers``, ``mesh_waterline_vertices``.

The variables ``diffraction_force``, ``excitation_force`` and ``Froude_Krylov_force`` are historical exceptions.
