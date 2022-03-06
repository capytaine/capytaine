
import logging
import capytaine as cpt

import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=160)
from scipy.linalg import block_diag

import meshmagick
import meshmagick.mesh as mm

# The hydrostatics module changed with version 3.0 of Meshmagick.
# This example uses the older version which is still available
# as 'hydrostatics_old' in recent versions of Meshmagick.
from packaging import version
if version.parse(meshmagick.__version__) < version.parse('3.0'):
    import meshmagick.hydrostatics as hs
else:
    import meshmagick.hydrostatics_old as hs


rho_water = 1025
g = 9.81

# Initialize floating body
r = 1.0
sphere = cpt.Sphere(
    radius=r,          # Dimension
    center=(0, 0, 0),    # Position
    nphi=20, ntheta=20,  # Fineness of the mesh
)
sphere.add_all_rigid_body_dofs()

# Visualize the body
# sphere.show()

# Create a hydrostatics body in Meshmagick
hsd = hs.Hydrostatics(mm.Mesh(sphere.mesh.vertices, sphere.mesh.faces),
					  cog=(0,0,0),
					  rho_water=rho_water,
					  grav=g).hs_data

# Inertial properties for neutrally buoyant constant density body
m = hsd['disp_mass']
I = np.array([[hsd['Ixx'], -1*hsd['Ixy'], -1*hsd['Ixz']],
              [-1*hsd['Ixy'], hsd['Iyy'], -1*hsd['Iyz']],
              [-1*hsd['Ixz'], -1*hsd['Iyz'], hsd['Izz']]])
M = block_diag(m, m, m, I)
sphere.mass = sphere.add_dofs_labels_to_matrix(M)
print(sphere.mass)

assert np.isclose(m, 1/2 * 4/3 * r**3 * np.pi * rho_water, rtol=1e-1)

# Hydrostatics
kHS = block_diag(0,0,hsd['stiffness_matrix'],0)
sphere.hydrostatic_stiffness = sphere.add_dofs_labels_to_matrix(kHS)

print(sphere.hydrostatic_stiffness)
assert np.isclose(kHS[2,2], r**2 * np.pi * rho_water * g, rtol=1e-1)
