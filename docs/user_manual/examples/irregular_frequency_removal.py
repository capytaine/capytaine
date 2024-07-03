import numpy as np
import capytaine as cpt

cpt.set_logging('INFO')

radius = 12.5
draught = 37.5
n_radius = 15
n_theta = 25
n_z = 30

# Define the range of frequencies to solve
omega_range = np.linspace(1.5, 2.25, 50)

solver = cpt.BEMSolver()

def compute_excitation_force(body, omega_range):
    problems = [
        cpt.DiffractionProblem(body=body, omega=omega)
        for omega in omega_range
    ]
    results = solver.solve_all(problems)
    data = cpt.assemble_dataset(results)
    return data['excitation_force']


# Initialize floating body by generating a geometric mesh
cylinder_mesh = cpt.mesh_vertical_cylinder(
    length=draught*2, radius=radius, center=(0, 0, 0), faces_max_radius=1.0,
    ).immersed_part()

### Body Without Lid (Body A)
body_1 = cpt.FloatingBody(mesh=cylinder_mesh)
body_1.add_translation_dof(name="Surge")

### Body With Lid at z = -3.699 (Body B)
lid_mesh = cylinder_mesh.generate_lid(z=-3.699)
body_2 = cpt.FloatingBody(mesh=cylinder_mesh, lid_mesh=lid_mesh)
body_2.add_translation_dof(name="Surge")

### Body With Lid at z = auto (Body C)
lid_mesh = cylinder_mesh.generate_lid(z=cylinder_mesh.lowest_lid_position(omega_max=omega_range.max()))
body_3 = cpt.FloatingBody(mesh=cylinder_mesh, lid_mesh=lid_mesh)
body_3.add_translation_dof(name="Surge")

force_1 = compute_excitation_force(body_1, omega_range)
force_2 = compute_excitation_force(body_2, omega_range)
force_3 = compute_excitation_force(body_3, omega_range)

# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
plt.plot(
    omega_range,
    np.abs(force_1.sel(influenced_dof='Surge')),
    marker='+',
    label='no lid'
)
plt.plot(
    omega_range,
    np.abs(force_2.sel(influenced_dof='Surge')),
    marker='*',
    label=r'$z = -3.699$m'
)
plt.plot(
    omega_range,
    np.abs(force_3.sel(influenced_dof='Surge')),
    marker='x',
    label='auto'
)

plt.xlabel(r'$\omega$')
plt.ylabel(r'$F_1 (abs)$')
plt.legend()
plt.tight_layout()
plt.show()
