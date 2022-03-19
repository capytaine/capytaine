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

    elongate_in_z_hs = body.get_hydrostatic_stiffnessij("elongate_in_z", "elongate_in_z", 
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


Output is
::

    The rotation dof Roll has been initialized around the origin of the domain (0, 0, 0).
    The rotation dof Pitch has been initialized around the origin of the domain (0, 0, 0).
    The rotation dof Yaw has been initialized around the origin of the domain (0, 0, 0).
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    Clipping slice_of_top_side_of_cylinder_375_mesh by Plane(normal=[0. 0. 1.], point=[0. 0. 0.]): all vertices are removed.
    grav:
        Capytaine  - 9.80665
        Meshmagick - 9.80665
    rho_water:
        Capytaine  - 1000
        Meshmagick - 1000
    cog:
        Capytaine  - [0, 0, 0]
        Meshmagick - [0 0 0]
    wet_surface_area:
        Capytaine  - 9.119266148961312
        Meshmagick - 9.119266148961316
    disp_volume:
        Capytaine  - 2.938926261462367
        Meshmagick - 2.9389262614623664
    disp_mass:
        Capytaine  - 2938.9262614623667
        Meshmagick - 2938.9262614623663
    buoyancy_center:
        Capytaine  - [ 4.12178644e-17 -1.96704181e-17 -5.00000000e-01]
        Meshmagick - [-2.36103028e-17 -3.77764845e-17 -5.00000000e-01]
    waterplane_center:
        Capytaine  - [1.77077271e-17 3.77764845e-17 0.00000000e+00]
        Meshmagick - [-6.29608076e-17 -3.14804038e-17  0.00000000e+00]
    waterplane_area:
        Capytaine  - 2.938926261462366
        Meshmagick - 2.9389262614623655
    transversal_metacentric_radius:
        Capytaine  - 0.22575292568380928
        Meshmagick - 0.23408474953124558
    longitudinal_metacentric_radius:
        Capytaine  - 0.22575292568380928
        Meshmagick - 0.23408474953124558
    transversal_metacentric_height:
        Capytaine  - -0.27424707431619055
        Meshmagick - -0.2659152504687545
    longitudinal_metacentric_height:
        Capytaine  - -0.27424707431619055
        Meshmagick - -0.2659152504687545
    stiffness_matrix:
        Capytaine  - [[ 2.88210212e+04  1.08875686e-12  5.44378431e-13]
    [ 1.08875686e-12 -7.90408075e+03 -6.80473039e-13]
    [ 5.44378431e-13 -6.80473039e-13 -7.90408075e+03]]
        Meshmagick - [[28821.02122197     0.             0.        ]
    [    0.         -7663.94907701     0.        ]
    [    0.             0.         -7663.94907701]]
    draught:
        Capytaine  - 1.0000000000000002
        Meshmagick - 1.0000000000000002
    length_at_waterline:
        Capytaine  - 2.0
        Meshmagick - 2.0
    breadth_at_waterline:
        Capytaine  - 1.9021130325903073
        Meshmagick - 1.9021130325903073
    length_overall_submerged:
        Capytaine  - 2.0
        Meshmagick - 2.0
    breadth_overall_submerged:
        Capytaine  - 1.9021130325903073
        Meshmagick - 1.9021130325903073
    inertia_matrix:
        Capytaine  - [[ 1.63731550e+03 -4.16333634e-14 -8.63303903e-15]
    [-4.16333634e-14  1.63731550e+03  1.11022302e-13]
    [-8.63303903e-15  1.11022302e-13  1.32840873e+03]]
        Meshmagick - [[ 1.66759990e+03 -1.56570115e-14 -2.77555756e-14]
    [-1.56570115e-14  1.66759990e+03 -5.55111512e-14]
    [-2.77555756e-14 -5.55111512e-14  1.37591564e+03]]


Verifying with Analytical Results
---------------------------------

Example code to verify with analytical results
::

    radius = 10
    cog = [0,0,0]
    body = cpt.Sphere(
        radius=radius,
        center=cog,
        nphi=300, ntheta=300,
    )
    body.add_all_rigid_body_dofs()
    density = 1000
    gravity = 9.80665

    
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

Output is 
::

    The rotation dof Roll has been initialized around the origin of the domain (0, 0, 0).
    The rotation dof Pitch has been initialized around the origin of the domain (0, 0, 0).
    The rotation dof Yaw has been initialized around the origin of the domain (0, 0, 0).
    wet_surface_area:
        Capytaine  - 628.2869509842606
        Meshmagick - 628.2869509842606
        Analytical - 628.3185307179587
    disp_volume:
        Capytaine  - 2094.18457402721
        Meshmagick - 2094.1845740271997
        Analytical - 2094.3951023931954
    buoyancy_center:
        Capytaine  - [-2.60089499e-16  7.21935327e-16 -3.74993146e+00]
        Meshmagick - [-2.17147694e-16 -1.52681972e-16 -3.74996573e+00]
        Analytical - [ 0.    0.   -3.75]
    waterplane_area:
        Capytaine  - 314.1362982503544
        Meshmagick - 314.13629825035423
        Analytical - 314.1592653589793
    transversal_metacentric_radius:
        Capytaine  - 3.7496573214020334
        Meshmagick - 3.749828657085805
        Analytical - 3.75
    longitudinal_metacentric_radius:
        Capytaine  - 3.749657321402033
        Meshmagick - 3.749828657085804
        Analytical - 3.75
    transversal_metacentric_height:
        Capytaine  - -0.0002741387393978556
        Meshmagick - -0.00013707282811603605
        Analytical - 0.0
    longitudinal_metacentric_height:
        Capytaine  - -0.0002741387393982997
        Meshmagick - -0.00013707282811692423
        Analytical - 0.0
    stiffness_matrix:
        Capytaine  - [[ 3.08062473e+06 -8.36165270e-10  1.11488703e-09]
    [-8.36165270e-10 -5.62996951e+03 -2.22977405e-09]
    [ 1.11488703e-09 -2.22977405e-09 -5.62996951e+03]]
        Meshmagick - [[ 3.08062473e+06  0.00000000e+00  0.00000000e+00]
    [ 0.00000000e+00 -2.81505578e+03  0.00000000e+00]
    [ 0.00000000e+00  0.00000000e+00 -2.81505578e+03]]
        Analytical - [[3080849.95963263       0.               0.        ]
    [      0.               0.               0.        ]
    [      0.               0.               0.        ]]