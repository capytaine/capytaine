============
Hydrostatics
============

Capytaine can compute some of the hydrostatic parameters of a given :code:`FloatingBody`. 


Integration
-----------

Here the integration approximated as the summation of the data function at face center multipled with respective face area

.. math::

    \iint_S f(x,y,z) dS \approx \sum_i^N f(x_i, y_i, z_i) \Delta S_i

where :math:`i` is the face index, :math:`(x_i, y_i, z_i)` is ith face center, and :math:`S_i` is ith face area. 

Hydrostatic Parameters
----------------------

For each hydrostatic parameter a separate method is available in Capytaine.

.. note::
    Before computing individual hydrostatic parameters, make sure to crop the body to only keep immersed.

::

    import capytaine as cpt
    import numpy as np

    cog = (0,0,-2)
    sphere = cpt.Sphere(radius=1.0, center=cog, name="my buoy")
    
    sphere.keep_immersed_part()

    hydrostatics = {}
    
    # :code:`get_volume()` returns volumes computes using x, y, z coordinates. 
    # This is similar to VOLX, VOLY, VOLZ in WAMIT.::
    hydrostatics["total_volumes"] = sphere.get_volumes() # [VOLX, VOLY, VOLZ]

    # Mean of the [VOLX, VOLY, VOLZ]
    hydrostatics["total_volume"] = sphere.get_volume()

    # Center of buoyancy
    hydrostatics["buoyancy_center"] = sphere.get_buoyancy_center()
    
    # Wet Surface Area
    hydrostatics["wet_surface_area"] = sphere.get_wet_surface_area()
    
    # Displaced Volume
    hydrostatics["disp_volume"] = sphere.get_volume()
    
    # Displaced Mass
    hydrostatics["disp_mass"] = sphere.get_mass()
    
    # Water Plane Center. Returns (0,0,0) for fully submerged bodies
    hydrostatics["waterplane_center"] = sphere.get_waterplane_center()
    
    # Water Plane Area. Returns 0.0 for fully submerged bodies
    hydrostatics["waterplane_area"] = sphere.get_waterplane_area()
    
    # Metacentric Parameters
    hydrostatics["transversal_metacentric_radius"] = sphere.get_bmt()
    hydrostatics["longitudinal_metacentric_radius"] = sphere.get_bml()
    hydrostatics["transversal_metacentric_height"] = sphere.get_gmt(cog=cog)
    hydrostatics["longitudinal_metacentric_height"] = sphere.get_gml(cog=cog)
    

Hydrostatic Stiffness
---------------------

The equation to compute the hydrostatic stiffness of a floating body is

.. math::

    C_{ij} = \iint_S (\hat{n} \cdot V_j) (w_i + z D_i)  dS
        
where :math:`\hat{n}` is surface normal, 

:math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and

:math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.


:code:`get_hydrostatic_stiffness` computes the hydrostatic stiffness using above equation directly from the DOFs (even for the rigid DOFs) and returns a (DOF count x DOF count) 2D matrix. ::  

    sphere.add_all_rigid_body_dofs()
    hydrostatics["stiffness_matrix"] = sphere.get_hydrostatic_stiffness()

    print(f"DOF count = {len(sphere.dofs)}")
    # DOF count = 6
    
    print(hydrostatics['stiffness_matrix'].shape)
    # (6, 6)


.. note::
    This method computes the hydrostatic stiffness assuming zero divergence. :math:`D_{i} = 0`. If :math:`D_i \neq 0`, input the divergence interpolated to face centers. 

::
  
    body.dofs["elongate_in_z"] = np.zeros_like(faces_centers)
    body.dofs["elongate_in_z"][:,2] = body.mesh.faces_centers[:,2]

    elongate_in_z_divergence = np.ones(body.mesh.faces_centers.shape[0])

    density = 1000
    gravity = 9.80665

    elongate_in_z_hs = body.get_hydrostatic_stiffnessij(body.dofs["elongate_in_z"], 
                                        body.dofs["elongate_in_z"], 
                                        divergence_i=elongate_in_z_divergence,
                                        density=density, gravity=gravity)

    analytical_hs = - density * gravity * (4 * body.volume() * body.get_buoyancy_center()[2])

    print( np.isclose(elongate_in_z_hs, analytical_hs) )
    # True


Interia Matrix
--------------

:code:`get_rigid_dof_mass` method computes 6 x 6 interia mass matrix of 6 rigid dofs. ::

    mass_matrix = body.get_rigid_dof_mass()

.. note::
    Unlike :code:`get_hydrostatic_stiffness`, the :code:`get_rigid_dof_mass` can only compute for 6 x 6 rigid interia mass. 

Compute all Hydrostatics
------------------------

Instead of computing each hydrostatic parameters, :code:`compute_hydrostatics` method computes all hydrostatic parameters and returns hydrostatic parameters :code:`dict`. 

.. note::
    No need to apply :code:`keep_immersed_part` to use :code:`compute_hydrostatics`.
    
::

    hydrostatics = body.compute_hydrostatics()

    print(hydrostatics.keys())

    # dict_keys(['grav', 'rho_water', 'cog', 'total_volume', 
    # 'total_volume_center', 'wet_surface_area', 'disp_volume', 
    # 'disp_mass', 'buoyancy_center', 'waterplane_center', 
    # 'waterplane_area', 'transversal_metacentric_radius', 
    # 'longitudinal_metacentric_radius', 'transversal_metacentric_height', 
    # 'longitudinal_metacentric_height', 'stiffness_matrix', 
    # 'length_overall', 'breadth_overall', 'depth', 'draught', 
    # 'length_at_waterline', 'breadth_at_waterline', 
    # 'length_overall_submerged', 'breadth_overall_submerged', 
    # 'inertia_matrix'])


Verifying with meshmagick library
---------------------------------

You can verify the results with meshmagick results
::

    import capytaine as cpt
    import numpy as np

    cog = [0,0,0]
    body = cpt.VerticalCylinder(
        length=2, radius=1.0,  # Dimensions
        center=cog,        # Position
        nr=10, nx=10, ntheta=10,   # Fineness of the mesh
    )
    body.add_all_rigid_body_dofs()

    capy_hsdb = body.compute_hydrostatics(cog=cog)
    capy_hsdb["stiffness_matrix"] = capy_hsdb["stiffness_matrix"][2:5,2:5]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"][3:,3:]


    import meshmagick.mesh as mmm
    import meshmagick.hydrostatics as mmhs

    body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)

    mm_hsdb = mmhs.compute_hydrostatics(body_mesh, np.array(cog), density, gravity)

    mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(rho_medium=density).inertia_matrix
    mm_hsdb["mesh"] = ""

    for var in capy_hsdb:
        if var in mm_hsdb:
            print(f"{var}:")
            print(f"    Capytaine  - {capy_hsdb[var]}")
            print(f"    Meshmagick - {mm_hsdb[var]}")




Verifying with Analytical Results
---------------------------------

Example code to verify with analytical results
::

    radius = 10

    body = cpt.Sphere(
        radius=radius,
        center=(0,0,0),
        nphi=50, ntheta=50,
    )
    body.add_all_rigid_body_dofs()
    # body.show()
    self=body
    density = 1000
    gravity = 9.80665

    cog = [0,0,0]
    capy_hsdb = body.compute_hydrostatics(density=density, gravity=gravity, cog=cog)
    capy_hsdb["stiffness_matrix"] = capy_hsdb["stiffness_matrix"][2:5,2:5]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"][3:,3:]


    import meshmagick.mesh as mmm
    import meshmagick.hydrostatics as mmhs

    # body.keep_immersed_part()
    body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)

    mm_hsdb = mmhs.compute_hydrostatics(body_mesh, np.array(cog), density, gravity)

    mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(rho_medium=density).inertia_matrix
    mm_hsdb["mesh"] = ""


    analytical = {}
    analytical["waterplane_area"] = np.pi*radius**2
    analytical["wet_surface_area"] = 2*np.pi*radius**2
    analytical["disp_volume"] = (2/3)*np.pi*radius**3
    analytical["interia_xx"] = np.pi*radius**4/4
    analytical["interia_yy"] = np.pi*radius**4/4
    analytical["interia_zz"] = np.pi*radius**4/2
    analytical["buoyancy_center"] = np.array([0,0,-analytical["interia_zz"] / (2*analytical["disp_volume"])])
    analytical["buoyancy_center"] = np.array([0,0,-3*radius/8])
    analytical["transversal_metacentric_radius"] = analytical["interia_xx"] / analytical["disp_volume"]
    analytical["longitudinal_metacentric_radius"] = analytical["interia_yy"] / analytical["disp_volume"]
    analytical["transversal_metacentric_height"] = analytical["transversal_metacentric_radius"] + analytical["buoyancy_center"][2] - cog[2]
    analytical["longitudinal_metacentric_height"] = analytical["longitudinal_metacentric_radius"] + analytical["buoyancy_center"][2] - cog[2]
    analytical["stiffness_matrix"] = density * gravity * np.array([
        [analytical["waterplane_area"], 0, 0],
        [0, analytical["disp_volume"] * analytical["transversal_metacentric_height"], 0],
        [0, 0, analytical["disp_volume"] * analytical["transversal_metacentric_height"]],
        ])

    for var in capy_hsdb:
        if var in analytical:
            print(f"{var}:")
            print(f"    Capytaine  - {capy_hsdb[var]}")
            print(f"    Meshmagick - {mm_hsdb[var]}")
            print(f"    Analytical - {analytical[var]}")
