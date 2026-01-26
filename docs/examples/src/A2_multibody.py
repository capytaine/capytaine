import numpy as np
import xarray as xr
import capytaine as cpt

cpt.set_logging('DEBUG')

# Define the first body
mesh_1 = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), resolution=(20, 20))
sphere = cpt.FloatingBody(
    mesh=mesh_1,
    lid_mesh=mesh_1.generate_lid(),
    dofs=cpt.rigid_body_dofs(only=['Surge', 'Heave'], rotation_center=(0, 0, 0)),
    center_of_mass=(0, 0, 0),
    name="sphere_1",
)

# Define the second body
mesh_2 = cpt.mesh_sphere(radius=0.5, center=(-2, -3, 0), resolution=(20, 20))
other_sphere = cpt.FloatingBody(
    mesh=mesh_2,
    lid_mesh=mesh_2.generate_lid(),
    dofs=cpt.rigid_body_dofs(only=['Surge', 'Heave'], rotation_center=(-2, -3, 0)),
    center_of_mass=(-2, -3, 0),
    name="sphere_2",
)

# Define the third body
mesh_3 = cpt.mesh_horizontal_cylinder(
    length=5.0, radius=1.0, center=(0.5, 3.0, -0.5), resolution=(3, 20, 20)
)
cylinder = cpt.FloatingBody(
    mesh=mesh_3,
    lid_mesh=mesh_3.generate_lid(),
    dofs=cpt.rigid_body_dofs(only=['Surge', 'Heave'], rotation_center=(0.5, 3.0, -0.5)),
    center_of_mass=(0.5, 3.0, -0.5),
    name="cylinder",
)

# Combine the three individual bodies into a single body.
all_bodies = cpt.Multibody([sphere, other_sphere, cylinder])

# all_bodies.mesh.immersed_part().show()

print("Merged body name:", all_bodies.name)
print("Merged body dofs:", list(all_bodies.dofs.keys()))
print("Multibody centers of mass:", all_bodies.center_of_mass)

# Compute multibody hydrostatics
M = all_bodies.compute_rigid_body_inertia()
K = all_bodies.compute_hydrostatic_stiffness()

with np.printoptions(precision=2, suppress=True):
    print(M)
    print(K)

# Compute multibody hydrodynamics
test_matrix = xr.Dataset(
    {
        "omega": np.linspace(0.5, 2.0, 2),
        "wave_direction": [0.0],
        "radiating_dof": list(all_bodies.dofs),
        "water_depth": [np.inf],
    }
)
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, all_bodies.immersed_part())

print(dataset)
