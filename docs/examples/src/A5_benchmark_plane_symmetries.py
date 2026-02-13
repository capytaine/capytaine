"""
This example benchmarks the performance of the mesh plane symmetry.
The same problems are solved, first on the full reference mesh, then using one symmetry plane, then using two symmetry planes.
"""

import numpy as np
import xarray as xr
import capytaine as cpt

cpt.set_logging("WARNING")

reference_mesh = cpt.mesh_parallelepiped(resolution=(30, 30, 30), name="reference_mesh")
print("Number of faces in full mesh:", reference_mesh.nb_faces)

# Cut the mesh to build symmetric mesh object
half_mesh = reference_mesh.clipped(origin=(0, 0, 0), normal=(1, 0, 0))
simple_symmetric_mesh = cpt.ReflectionSymmetricMesh(
    half=half_mesh, plane="yOz", name="simple_symmetric_mesh"
)

quarter_mesh = half_mesh.clipped(origin=(0, 0, 0), normal=(0, 1, 0))
nested_symmetric_mesh = cpt.ReflectionSymmetricMesh(
    half=cpt.ReflectionSymmetricMesh(half=quarter_mesh, plane="xOz"),
    plane="yOz",
    name="nested_symmetric_mesh",
)

for mesh in [reference_mesh, simple_symmetric_mesh, nested_symmetric_mesh]:
    body = cpt.FloatingBody(
        mesh=mesh,
        lid_mesh=mesh.generate_lid(),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.25)),
        name=mesh.name,
    )

    test_matrix = xr.Dataset(
        {
            "omega": np.linspace(0.5, 2.0, 5),
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(body.dofs),
        }
    )

    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(
        test_matrix,
        body.immersed_part(),
        n_threads=1,  # No parallel resolution for cleaner benchmark
    )
    print()
    print("--->", body.name)
    print(
        "Added mass",
        dataset.added_mass.sel(radiating_dof="Heave", influenced_dof="Heave").values,
    )
    print("Computation time", solver.timer_summary())
    print()
