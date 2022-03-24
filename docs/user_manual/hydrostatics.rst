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

    sphere.center_of_mass = cog

    sphere.keep_immersed_part()

    hydrostatics = {}

    # :code:`volumes` returns volumes computes using x, y, z coordinates. 
    # This is similar to VOLX, VOLY, VOLZ in WAMIT.::
    hydrostatics["total_volumes"] = sphere.volumes # [VOLX, VOLY, VOLZ]

    # Center of buoyancy
    hydrostatics["buoyancy_center"] = sphere.center_of_buoyancy
    # Wet Surface Area
    hydrostatics["wet_surface_area"] = sphere.wet_surface_area

    # Displaced Volume
    hydrostatics["disp_volume"] = sphere.volume

    # Displaced Mass
    hydrostatics["disp_mass"] = sphere.disp_mass()

    # Water Plane Center. Returns (0,0,0) for fully submerged bodies
    hydrostatics["waterplane_center"] = sphere.waterplane_center

    # Water Plane Area. Returns 0.0 for fully submerged bodies
    hydrostatics["waterplane_area"] = sphere.waterplane_area

    # Metacentric Parameters
    hydrostatics["transversal_metacentric_radius"] = sphere.bmt
    hydrostatics["longitudinal_metacentric_radius"] = sphere.bml
    hydrostatics["transversal_metacentric_height"] = sphere.gmt
    hydrostatics["longitudinal_metacentric_height"] = sphere.gml
    

Hydrostatic Stiffness
---------------------

The equation to compute the hydrostatic stiffness of a floating body is

.. math::

    C_{ij} = \iint_S (\hat{n} \cdot V_j) (w_i + z D_i)  dS
        
where :math:`\hat{n}` is surface normal, 

:math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and

:math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.


:code:`hydrostatic_stiffness_xr` computes the hydrostatic stiffness using above equation directly from the DOFs (even for the rigid DOFs) and returns a (DOF count x DOF count) 2D matrix. ::  

    sphere.add_all_rigid_body_dofs()
    hydrostatics["stiffness_matrix"] = sphere.hydrostatic_stiffness_xr()

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

    elongate_in_z_hs = body.each_hydrostatic_stiffness("elongate_in_z", "elongate_in_z", 
                                        divergence_i=elongate_in_z_divergence,
                                        density=density, gravity=gravity)

    analytical_hs = - density * gravity * (4 * body.volume * body.center_of_buoyancy[2])

    print( np.isclose(elongate_in_z_hs, analytical_hs) )
    # True


Interia Matrix
--------------

:code:`rigid_dof_mass` method computes 6 x 6 interia mass matrix of 6 rigid dofs. ::

    mass_matrix = body.rigid_dof_mass()

.. note::
    Unlike :code:`hydrostatic_stiffness_xr`, the :code:`rigid_dof_mass` can only compute for 6 x 6 rigid interia mass. 

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


Verifying with Meshmagick and Analytical Results
------------------------------------------------

Example code to verify with Meshmagick and Analytical results
::

    import capytaine as cpt
    import numpy as np
    import meshmagick.mesh as mmm
    import meshmagick.hydrostatics as mmhs

    radius = 10
    cog = (0,0,0)
    body = cpt.Sphere(
        radius=radius,
        center=cog,
        nphi=100, ntheta=100,
    )
    body.center_of_mass = cog

    body.keep_immersed_part()
    body.add_all_rigid_body_dofs()
    # body.show()
    self=body

    density = 1000
    gravity = 9.80665

    capy_hsdb = body.compute_hydrostatics(density=density, gravity=gravity)

    stiff_compare_dofs = ["Heave", "Roll", "Pitch"]
    capy_hsdb["stiffness_matrix"] = capy_hsdb["stiffness_matrix"].sel(
        influenced_dof=stiff_compare_dofs, radiating_dof=stiff_compare_dofs
        ).values

    mass_compare_dofs = ["Roll", "Pitch", "Yaw"]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"].sel(
        influenced_dof=mass_compare_dofs, radiating_dof=mass_compare_dofs
        ).values


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
  
    wet_surface_area:
        Capytaine  - 628.0343659038494
        Meshmagick - 628.0343659038496
        Analytical - 628.3185307179587
    disp_volume:
        Capytaine  - 2092.5009287939088
        Meshmagick - 2092.5009287939115
        Analytical - 2094.3951023931954
    waterplane_area:
        Capytaine  - 313.95259764656686
        Meshmagick - 313.95259764656674
        Analytical - 314.1592653589793
    transversal_metacentric_radius:
        Capytaine  - 3.7469169327091647
        Meshmagick - 3.748458229464248
        Analytical - 3.75
    longitudinal_metacentric_radius:
        Capytaine  - 3.7469169327091643
        Meshmagick - 3.748458229464248
        Analytical - 3.75
    transversal_metacentric_height:
        Capytaine  - -0.002466140909095582
        Meshmagick - -0.0012332946572213288
        Analytical - 0.0
    longitudinal_metacentric_height:
        Capytaine  - -0.0024661409090960262
        Meshmagick - -0.0012332946572213288
        Analytical - 0.0
    stiffness_matrix:
        Capytaine  - [[ 3.07882324e+06 -1.11488703e-09  0.00000000e+00]
     [-1.11488703e-09 -5.06062577e+04  2.22977405e-09]
     [ 0.00000000e+00  2.22977405e-09 -5.06062577e+04]]
        Meshmagick - [[3078823.2417107        0.               0.        ]
     [      0.          -25307.72957091       0.        ]
     [      0.               0.          -25307.72957091]]
        Analytical - [[3080849.95963263       0.               0.        ]
     [      0.               0.               0.        ]
     [      0.               0.               0.        ]]