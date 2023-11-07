import numpy as np

import capytaine as cpt
import meshmagick.mesh
import meshmagick.hydrostatics

radius = 10
cog = np.array((0, 0, 0))
mesh = cpt.mesh_sphere(radius=radius, center=cog, resolution=(100, 100))
body = cpt.FloatingBody(mesh=mesh, center_of_mass=cog)

body.add_all_rigid_body_dofs()
body = body.immersed_part()

density = 1000
gravity = 9.81

capy_hsdb = body.compute_hydrostatics(rho=density, g=gravity)

stiff_compare_dofs = ["Heave", "Roll", "Pitch"]
capy_hsdb["stiffness_matrix"] = capy_hsdb["hydrostatic_stiffness"].sel(
    influenced_dof=stiff_compare_dofs, radiating_dof=stiff_compare_dofs
    ).values

mass_compare_dofs = ["Roll", "Pitch", "Yaw"]
capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"].sel(
    influenced_dof=mass_compare_dofs, radiating_dof=mass_compare_dofs
    ).values


mm_mesh = meshmagick.mesh.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)

mm_hsdb = meshmagick.hydrostatics.compute_hydrostatics(mm_mesh, cog=cog, rho_water=density, grav=gravity)

mm_hsdb["inertia_matrix"] = mm_mesh.eval_plain_mesh_inertias(rho_medium=density).inertia_matrix
mm_hsdb["center_of_buoyancy"] = mm_hsdb["buoyancy_center"]


analytical = {}
analytical["waterplane_area"] = np.pi*radius**2
analytical["wet_surface_area"] = 2*np.pi*radius**2
analytical["disp_volume"] = (2/3)*np.pi*radius**3
analytical["center_of_buoyancy"] = np.array([0,0,-3*radius/8])
inertia_xx = np.pi*radius**4/4
inertia_yy = np.pi*radius**4/4
inertia_zz = np.pi*radius**4/2
analytical["transversal_metacentric_radius"] = inertia_xx / analytical["disp_volume"]
analytical["longitudinal_metacentric_radius"] = inertia_yy / analytical["disp_volume"]
analytical["transversal_metacentric_height"] = analytical["transversal_metacentric_radius"] + analytical["center_of_buoyancy"][2] - cog[2]
analytical["longitudinal_metacentric_height"] = analytical["longitudinal_metacentric_radius"] + analytical["center_of_buoyancy"][2] - cog[2]

analytical["stiffness_matrix"] = density * gravity * np.array([
    [analytical["waterplane_area"], 0, 0],
    [0, analytical["disp_volume"] * analytical["transversal_metacentric_height"], 0],
    [0, 0, analytical["disp_volume"] * analytical["transversal_metacentric_height"]],
    ])

vars_to_be_displayed = capy_hsdb.keys() & (mm_hsdb.keys() | analytical.keys())
for var in vars_to_be_displayed:
    print(f"{var}:")
    print(f"    Capytaine  :  {capy_hsdb[var]}")
    if var in mm_hsdb:
        print(f"    Meshmagick :  {mm_hsdb[var]}")
    if var in analytical:
        print(f"    Analytical :  {analytical[var]}")
