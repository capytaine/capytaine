============
Hydrostatics
============

Capytaine can compute some of the hydrostatic parameters of a given :code:`FloatingBody`.


Hydrostatic parameters
----------------------

For each hydrostatic parameter a separate method is available in Capytaine.
Some of them may require the definition of the center of mass of the body.
It can be done by setting the attribute `center_of_mass` as in the example below.

::

    import capytaine as cpt
    import numpy as np

    rigid_sphere = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(radius=1.0, center=(0, 0, 0)),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.3)),
            center_of_mass=(0, 0, -0.3)
            )

    print("Volume:", rigid_sphere.volume)
    print("Center of buoyancy:", rigid_sphere.center_of_buoyancy)
    print("Wet surface area:", rigid_sphere.wet_surface_area)
    print("Displaced mass:", rigid_sphere.disp_mass(rho=1025))
    print("Waterplane center:", rigid_sphere.waterplane_center)
    print("Waterplane area:", rigid_sphere.waterplane_area)
    print("Metacentric parameters:",
        rigid_sphere.transversal_metacentric_radius,
        rigid_sphere.longitudinal_metacentric_radius,
        rigid_sphere.transversal_metacentric_height,
        rigid_sphere.longitudinal_metacentric_height)


Hydrostatic stiffness
---------------------

The method :meth:`~capytaine.bodies.FloatingBody.compute_hydrostatic_stiffness`
computes the hydrostatic stiffness and returns a (DOF count x DOF count) 2D
matrix as an :code:`xarray.DataArray`. ::

    hs = rigid_sphere.compute_hydrostatic_stiffness()

Capytaine has two built-in methods to compute the hydrostatic stiffness:

- for rigid-body degrees of freedom, exact analytical expressions are available.

- for generalized degrees of freedom (e.g. elastic bodies), the following expression is used

.. math::

    C_{ij} = \iint_S (\hat{n} \cdot V_j) (w_i + z D_i)  dS

where :math:`\hat{n}` is the surface normal,
:math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is the dof vector and
:math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.

Note that the latter method is an approximation and does not recover the former
expression for rigid-body dofs, even when passing a dof vector and a divergence
corresponding to the rigid-body dof.

.. warning::
   Currently, the detection of a rigid-body degree of freedom to use the exact
   analytical expression is not very robust. Changing the name of the dof or
   combining several rigid bodies usually makes Capytaine use the approximate
   method for generalize dofs instead.
   This is planned to be fixed in a later release of Capytaine.


If the divergence is not provided by the user, zero is used as an approximation.
Here is an example of setting the divergence of the dof::

    mesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0)).immersed_part()
    dof = np.array([(0, 0, z) for (x, y, z) in elastic_sphere.mesh.faces_centers])
    elastic_sphere = cpt.FloatingBody(
            mesh=mesh,
            center_of_mass=(0, 0, -0.3),
            dofs={"elongate_in_z": dof},
            )

    dofs_divergence = {"elongate_in_z": np.ones(elastic_sphere.mesh.nb_faces)}

    density = 1000.0
    gravity = 9.81
    elongate_in_z_hs = elastic_sphere.compute_hydrostatic_stiffness(divergence=dofs_divergence, rho=density, g=gravity)

    analytical_hs = - density * gravity * (4 * elastic_sphere.volume * elastic_sphere.center_of_buoyancy[2])
    print(np.isclose(elongate_in_z_hs.values[0, 0], analytical_hs))
    # True


.. warning::
   In the example above, all computations are done on a mesh clipped with
   ``immersed_part()``.
   Defining a divergence for generalized dofs is one of the only cases where
   the user should make sure that the body contains only immersed panels, and
   that the divergence field is defined only on these panels.
   It is planned to remove this requirement in a future release of Capytaine.

If the mass is not specified (as in the examples above), the body is assumed to
be in buoyancy equilibrium. Its mass is the mass of the displaced volume of
water.

Non-neutrally buoyant bodies are partially implemented (only for a single rigid body).
In this case, the mass of the body can be given by setting the :code:`mass`
attribute of the :code:`FloatingBody`::

    mesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0))
    rigid_sphere = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.3)),
            center_of_mass=(0, 0, -0.3),
            mass=0.8*mesh.disp_mass(rho=1000),
            )
    hs = rigid_sphere.compute_hydrostatic_stiffness()


Inertia matrix
--------------

The method :meth:`~capytaine.bodies.FloatingBody.compute_rigid_body_inertia` is
able to computes the 6 x 6 inertia matrix of a body with 6 rigid dofs.
The inertia coefficient of other degrees of freedom are filled with :code:`NaN` by default.

::

    M = rigid_sphere.compute_rigid_body_inertia()


As for the hydrostatic stiffness, the mass is assumed to be the displaced mass
of water, unless a :code:`mass` attribute has been specified.

A custom matrix can be provided. For consistency with the data computed with
Capytaine, it is recommended to wrap it in a :code:`xarray.DataArray` with dof
names as labels::

    elastic_sphere.inertia_matrix = elastic_sphere.add_dofs_labels_to_matrix(np.array([[1000.0]]))


Compute all hydrostatics parameters
-----------------------------------

Instead of computing each hydrostatic parameters individually, :code:`compute_hydrostatics` returns a :code:`dict` containing all hydrostatic parameters.

::

    hydrostatics = rigid_sphere.compute_hydrostatics()

    print(hydrostatics.keys())
    # dict_keys(['g', 'rho', 'center_of_mass', 'wet_surface_area', 'disp_volumes',
    # 'disp_volume', 'disp_mass', 'center_of_buoyancy', 'waterplane_center',
    # 'waterplane_area', 'transversal_metacentric_radius',
    # 'longitudinal_metacentric_radius' , 'transversal_metacentric_height',
    # 'longitudinal_metacentric_height', 'hydrostatic_stiffness',
    # 'length_overall', 'breadt h_overall', 'depth', 'draught',
    # 'length_at_waterline', 'breadth_at_waterline',
    # 'length_overall_submerged', 'breadth_overall_submerged', 'inertia_matrix'])
