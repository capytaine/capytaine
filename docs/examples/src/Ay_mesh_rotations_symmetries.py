import numpy as np
import xarray as xr
import capytaine as cpt
from capytaine.new_meshes.meshes import to_new_mesh
from capytaine.new_meshes.symmetric_meshes import ReflectionSymmetricMesh, RotationSymmetricMesh

# sphere = to_new_mesh(cpt.mesh_sphere(resolution=(20, 20), center=(2.0, 0.0, 0.0), name="sphere")).immersed_part()
sphere = to_new_mesh(cpt.mesh_parallelepiped(resolution=(10, 10, 10), center=(2.0, 0.0, 0.0))).immersed_part()
half_sphere = sphere.clipped(origin=(0, 0, 0), normal=(0.0, 1.0, 0.0))

symmetric_mesh = RotationSymmetricMesh(wedge=sphere, n=4)
nested_symmetric_mesh_1 = ReflectionSymmetricMesh(RotationSymmetricMesh(wedge=half_sphere, n=4), plane="xOz")
nested_symmetric_mesh_2 = RotationSymmetricMesh(wedge=ReflectionSymmetricMesh(half_sphere, plane="xOz"), n=4)

# symmetric_mesh.show()
# nested_symmetric_mesh_1.show()
# nested_symmetric_mesh_2.show()

for mesh in [
    symmetric_mesh.merged(),
    symmetric_mesh,
    nested_symmetric_mesh_1,
    nested_symmetric_mesh_2
]:
    body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.25)),
        name=mesh.name,
        )

    test_matrix = xr.Dataset({
            "omega": np.linspace(0.5, 2.0, 5),
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(body.dofs),
            "water_depth": [8.0],
            "rho": [1025],
            })

    solver = cpt.BEMSolver(method='direct')
    dataset = solver.fill_dataset(
        test_matrix,
        body.immersed_part(),
        n_threads=1,  # No parallel resolution for cleaner benchmark
        progress_bar=False,
    )

    print(mesh.nb_faces)
    print("Added mass", dataset.added_mass.sel(radiating_dof='Surge', influenced_dof='Surge').values)
    print("Diffraction", dataset.diffraction_force.real.sel(wave_direction=0.0, influenced_dof='Surge').values)
    print(solver.timer_summary())

    S = solver.engine.last_computed_matrices[0]
    print(type(S), end='')
    try:
        print([type(b) for b in S.blocks])
    except Exception:
        print()
    print()
