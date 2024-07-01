import numpy as np
import capytaine as cpt

radius  =12.5; draught = 37.5
nRadius = 10; nTheta = 25; nZ = 30

# Initialize floating body by generating a geometric mesh
cylinderMesh = cpt.mesh_vertical_cylinder(
    length=draught* 2, radius=radius,  # Dimensions
    center=(0, 0, 0),        # Position
    resolution=(nRadius, nTheta, nZ)    # Fineness of the mesh
    )

lid = cylinderMesh.generate_lid(z=-5)
body = cpt.FloatingBody(mesh=cylinderMesh, lid_mesh=lid, dofs=cpt.rigid_body_dofs())
body.keep_immersed_part(free_surface=0.0,water_depth=np.infty)

body.show()
