============
Hydrostatics
============

Capytaine can compute some of the hydrostatic parameters of a given :code:`FloatingBody`.


.. note::
    Here the integration of a function over the immersed hull is approximated as the summation of the data function at face centers multiplied by the respective face areas.

    .. math::

        \iint_S f(x,y,z) dS \approx \sum_i^N f(x_i, y_i, z_i) \Delta S_i

    where :math:`i` is the face index, :math:`(x_i, y_i, z_i)` is :math:`i`th face center, and :math:`S_i` is :math:`i`th face area.

Hydrostatic parameters
----------------------

For each hydrostatic parameter a separate method is available in Capytaine.
Some of them may require the definition of the center of mass of the body.
It can be done by setting the attribute `center_of_mass` as in the example below.

.. note::
    Before computing individual hydrostatic parameters, make sure to crop the body to only keep the immersed part.

::

    import capytaine as cpt
    import numpy as np

    rigid_sphere = cpt.Sphere(radius=1.0, center=(0, 0, 0), name="my buoy")
    rigid_sphere.center_of_mass = (0, 0, -0.3)
    rigid_sphere.add_all_rigid_body_dofs()
    rigid_sphere.keep_immersed_part()

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

The equation to compute the hydrostatic stiffness of a floating body is

.. math::

    C_{ij} = \iint_S (\hat{n} \cdot V_j) (w_i + z D_i)  dS

where :math:`\hat{n}` is the surface normal,

:math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is the dof vector and

:math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.

The method :meth:`~capytaine.bodies.FloatingBody.compute_hydrostatic_stiffness`
computes the hydrostatic stiffness and returns a (DOF count x DOF count) 2D
matrix as an :code:`xarray.DataArray`. ::

    rigid_sphere.hydrostatic_stiffness = rigid_sphere.compute_hydrostatic_stiffness()


.. note::
   For rigid body dofs, the exact formula above can be evaluated.
   However, in general, the divergence :math:`D_i` of an arbitrary dofs is not known.
   User can pass a value for `D_i` as an input to :code:`compute_hydrostatics`.
   Otherwise the code assumes zero divergence :math:`D_{i} = 0`.

::

    elastic_sphere = cpt.Sphere(radius=1.0, center=(0, 0, 0), name="my buoy")
    elastic_sphere.center_of_mass = (0, 0, -0.3)
    elastic_sphere.keep_immersed_part()

    elastic_sphere.dofs["elongate_in_z"] = np.array([(0, 0, z) for (x, y, z) in elastic_sphere.mesh.faces_centers])

    dofs_divergence = {"elongate_in_z": np.ones(elastic_sphere.mesh.nb_faces)}

    density = 1000.0
    gravity = 9.81
    elongate_in_z_hs = elastic_sphere.compute_hydrostatic_stiffness(divergence=dofs_divergence, rho=density, g=gravity)

    analytical_hs = - density * gravity * (4 * elastic_sphere.volume * elastic_sphere.center_of_buoyancy[2])

    print(np.isclose(elongate_in_z_hs.values[0, 0], analytical_hs))
    # True


If the mass is not specified (as in the examples above), the body is assumed to
be in buoyancy equilibrium. It's mass is the mass of the displaced volume of
water.

Non-neutrally buoyant bodies are partially implemented (only for a single rigid body).
In this case, the mass of the body can be given by setting the :code:`mass`
attribute of the :code:`FloatingBody`::

    rigid_sphere.mass = 0.8 * rigid_sphere.disp_mass(rho=1000)
    rigid_sphere.hydrostatic_stiffness = rigid_sphere.compute_hydrostatic_stiffness()


Inertia matrix
--------------

The method :meth:`~capytaine.bodies.FloatingBody.compute_rigid_body_inertia` is
able to computes the 6 x 6 inertia matrix of a body with 6 rigid dofs.
The inertia coefficient of other degrees of freedom are filled with :code:`NaN` by default.

::

    rigid_sphere.inertia_matrix = elastic_sphere.compute_rigid_body_inertia()


As for the hydrostatic stiffness, the mass is assumed to be the displaced mass
of water, unless a :code:`mass` attribute has been specified.

A custom matrix can be provided. For consistency with the data computed with
Capytaine, it is recommended to wrap it in a :code:`xarray.DataArray` with dof
names as labels::

    elastic_sphere.inertia_matrix = elastic_sphere.add_dofs_labels_to_matrix(np.array([[1000.0]]))


Compute all hydrostatics parameters
-----------------------------------

Instead of computing each hydrostatic parameters individually, :code:`compute_hydrostatics` returns a :code:`dict` containing all hydrostatic parameters.

.. note::
    No need to apply :code:`keep_immersed_part` to use :code:`compute_hydrostatics`.

::

    hydrostatics = rigid_sphere.compute_hydrostatics()

    print(hydrostatics.keys())
    # dict_keys(['g', 'rho', 'center_of_mass', 'wet_surface_area', 'disp_volumes',
    # 'disp_volume', 'disp_mass', 'center_ of_buoyancy', 'waterplane_center',
    # 'waterplane_area', 'transversal_metacentric_radius',
    # 'longitudinal_metacentric_radius' , 'transversal_metacentric_height',
    # 'longitudinal_metacentric_height', 'hydrostatic_stiffness',
    # 'length_overall', 'breadt h_overall', 'depth', 'draught',
    # 'length_at_waterline', 'breadth_at_waterline',
    # 'length_overall_submerged', 'breadth_overall_submerged', 'inertia_matrix'])

