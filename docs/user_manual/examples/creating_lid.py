#!/usr/bin/env python

import numpy as np
import capytaine as cpt

# cpt.set_logging('INFO')
######### GENERATE CYLINDER AND LID  ###################
'''
The vertical cylinder properties: 
radius  = 12.5m; draught = 37.5m; depth = infinity
nRadius = 40; nTheta = 30; nz = 20
location of center at x,y = (0,15) m
'''
radius  =12.5; draught = 37.5
nRadius = 15; nTheta = 25; nZ = 30

# Initialize floating body by generating a geometric mesh
cylinderMesh = cpt.mesh_vertical_cylinder(
    length=draught* 2, radius=radius,  # Dimensions
    center=(0, 0, 0),        # Position
    resolution=(nRadius, nTheta, nZ)    # Fineness of the mesh
    )
cylinder = cpt.FloatingBody(cylinderMesh)
cylinder.keep_immersed_part(free_surface=0.0,water_depth=np.infty)
lid = cylinder.generate_lid(z=-5,info=False)

body = cpt.FloatingBody(mesh=cylinder.mesh, lid_mesh=lid, dofs=cpt.rigid_body_dofs())
body.show()

