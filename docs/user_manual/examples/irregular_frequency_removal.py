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
based on On the Irregular Frequencies Appearing in Wave Diffraction-Radiation Solutions 
Malenica and Chen 1998 - OMAE
'''
radius  =12.5; draught = 37.5
nRadius = 15; nTheta = 25; nZ = 30
# Define the range of frequencies to solve 
omega_range = np.linspace(0.5, 2.25, 300)

def solve_excitation_basic(body,omega_range):
    # Calling solver 
    solver = cpt.BEMSolver()
    
    # Solution without lid
    problems = [
        cpt.DiffractionProblem(body=body,
                            omega=omega)
        for omega in omega_range
    ]
    results = solver.solve_all(problems,progress_bar=True)
    data = cpt.assemble_dataset(results)
    excitation_force = data['diffraction_force'] + data['Froude_Krylov_force']
    return excitation_force



# Initialize floating body by generating a geometric mesh
cylinderMesh = cpt.mesh_vertical_cylinder(
    length=draught* 2, radius=radius,  # Dimensions
    center=(0, 0, 0),        # Position
    resolution=(nRadius, nTheta, nZ)    # Fineness of the mesh
    )

### Body Without Lid (Body A)
bodyA = cpt.FloatingBody(mesh=cylinderMesh)
bodyA.add_translation_dof(name="Surge")
bodyA.keep_immersed_part(free_surface=0.0,water_depth=np.infty)

### Body With Lid at z = -3.699 (Body B)
lid = cylinderMesh.generate_lid(z=-3.699,info=False)
bodyB = cpt.FloatingBody(mesh=cylinderMesh, 
                        lid_mesh=lid)
bodyB.add_translation_dof(name="Surge")
bodyB.keep_immersed_part(free_surface=0.0,water_depth=np.infty)


### Body With Lid at z = auto (Body C)
lid = cylinderMesh.generate_lid(z='auto',omega_max =omega_range[-1],info=False)
bodyC = cpt.FloatingBody(mesh=cylinderMesh, 
                        lid_mesh=lid)
bodyC.add_translation_dof(name="Surge")
bodyC.keep_immersed_part(free_surface=0.0,water_depth=np.infty)

forceA = solve_excitation_basic(bodyA,omega_range)
forceB = solve_excitation_basic(bodyB,omega_range)
forceC = solve_excitation_basic(bodyC,omega_range)

# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
plt.plot(
    omega_range,
    np.abs(forceA.sel(influenced_dof='Surge')),
    marker='+',label = r'no lid'
)
plt.plot(
    omega_range,
    np.abs(forceB.sel(influenced_dof='Surge')),
    marker='*',label = r'$z =-3.699$m'
)
plt.plot(
    omega_range,
    np.abs(forceC.sel(influenced_dof='Surge')),
    marker='x',label = r'$z =auto$'
)

plt.xlabel(r'$\omega$')
plt.ylabel(r'$F_1 (abs)$')
plt.legend()
plt.tight_layout()
plt.show()

