import numpy as np
import capytaine as cpt

cpt.set_logging()

# Generating a geometric mesh
mesh = cpt.mesh_horizontal_cylinder(
    length=10.0,
    radius=1.0,
    center=(0, 0, -2,),
    resolution=(10, 20, 30,)
)

# Define a rigid body using this mesh
cylinder = cpt.FloatingBody(
    mesh=mesh,
    dofs=cpt.rigid_body_dofs(),
    name="floating cylinder"
)

# Define the range of frequencies as a Numpy array
omega_range = np.linspace(0.1, 3.0, 10)

# Set up the problems: we will solve a radiation problem for each
# degree of freedom of the body and for each frequency in the
# frequency range.
problems = [
    cpt.RadiationProblem(body=cylinder, radiating_dof=dof, omega=omega)
    for dof in cylinder.dofs
    for omega in omega_range
]
# Water density, gravity and water depth have not been specified.
# Default values are used.

# Solve all radiation problems
solver = cpt.BEMSolver()
results = solver.solve_all(problems)

# Gather the computed added mass into a labelled array.
data = cpt.assemble_dataset(results)

# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
for dof in cylinder.dofs:
    plt.plot(
        omega_range,
        data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )
plt.xlabel('omega')
plt.ylabel('added mass')
plt.legend()
plt.tight_layout()
plt.show()
