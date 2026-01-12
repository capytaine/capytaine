import numpy as np
import xarray as xr
import capytaine as cpt

reference_mesh = cpt.mesh_parallelepiped(resolution=(40, 40, 40), name="reference_mesh")
print(reference_mesh.nb_faces)
half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
quarter_mesh = half_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(0, 1, 0)))

old_simple_symmetric_mesh = cpt.ReflectionSymmetricMesh(half=half_mesh, plane=cpt.yOz_Plane, name="old_simple_symmetric_mesh")
old_nested_symmetric_mesh = cpt.ReflectionSymmetricMesh(half=cpt.ReflectionSymmetricMesh(half=quarter_mesh, plane=cpt.xOz_Plane), plane=cpt.yOz_Plane, name="old_nested_symmetric_mesh")

from capytaine.new_meshes.meshes import to_new_mesh
from capytaine.new_meshes.symmetric_meshes import ReflectionSymmetricMesh
new_simple_symmetric_mesh = ReflectionSymmetricMesh(half=to_new_mesh(half_mesh), plane="yOz", name="new_simple_symmetric_mesh")
new_nested_symmetric_mesh = ReflectionSymmetricMesh(half=ReflectionSymmetricMesh(half=to_new_mesh(quarter_mesh), plane="xOz"), plane="yOz", name="new_nested_symmetric_mesh")

# new_simple_symmetric_mesh.show()

for mesh in [reference_mesh, old_simple_symmetric_mesh, old_nested_symmetric_mesh, new_simple_symmetric_mesh, new_nested_symmetric_mesh]:
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
