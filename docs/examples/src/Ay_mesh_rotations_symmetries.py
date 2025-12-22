import numpy as np
import xarray as xr
import capytaine as cpt
from capytaine.new_meshes.meshes import to_new_mesh
from capytaine.new_meshes.symmetric_meshes import ReflectionSymmetricMesh, RotationSymmetricMesh

reference_mesh = cpt.mesh_sphere(resolution=(60, 60), name="reference_mesh")
new_reference_mesh = to_new_mesh(reference_mesh)
wedge = new_reference_mesh.extract_wedge(n=5)
half_wedge = wedge.rotated_z(-np.pi/5).clipped(origin=(0, 0, 0), normal=(0.0, 1.0, 0.0))

print(reference_mesh.nb_faces)

symmetric_mesh = RotationSymmetricMesh(wedge=wedge, n=5)
nested_symmetric_mesh_1 = ReflectionSymmetricMesh(RotationSymmetricMesh(wedge=half_wedge, n=5), plane="xOz")

# Currently not fully implemented, use the equivalent solution above
nested_symmetric_mesh_2 = RotationSymmetricMesh(wedge=ReflectionSymmetricMesh(half_wedge, plane="xOz"), n=5)

# symmetric_mesh.show()
# nested_symmetric_mesh_1.show()
# nested_symmetric_mesh_2.show()

for mesh in [
    reference_mesh,
    symmetric_mesh,
    nested_symmetric_mesh_1,
    # nested_symmetric_mesh_2
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

    solver = cpt.BEMSolver(method='indirect')
    dataset = solver.fill_dataset(
        test_matrix,
        body.immersed_part(),
        n_threads=1,  # No parallel resolution for cleaner benchmark
    )
    print(mesh.name)
    print(dataset.added_mass.sel(radiating_dof='Heave', influenced_dof='Heave').values)
    print(solver.timer_summary())
