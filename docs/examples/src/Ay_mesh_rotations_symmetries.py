import numpy as np
import xarray as xr
import capytaine as cpt
from capytaine.new_meshes.meshes import to_new_mesh
from capytaine.new_meshes.symmetric_meshes import RotationSymmetricMesh

single_pillar = cpt.mesh_parallelepiped(resolution=(10, 10, 10), center=(2.0, 0.0, 0.0)).immersed_part()
single_lid_mesh = single_pillar.generate_lid()
symmetric_mesh = RotationSymmetricMesh(wedge=to_new_mesh(single_pillar), n=4)
symmetric_lid_mesh = RotationSymmetricMesh(wedge=to_new_mesh(single_lid_mesh), n=4)

symmetric_body = cpt.FloatingBody(
        mesh=symmetric_mesh,
        lid_mesh=symmetric_lid_mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        )

full_body = cpt.FloatingBody(
        mesh=symmetric_mesh.merged(),
        lid_mesh=symmetric_lid_mesh.merged(),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        )

for body in [symmetric_body, full_body]:
    test_matrix = xr.Dataset({
            "omega": np.linspace(0.5, 2.0, 5),
            "wave_direction": np.linspace(0, np.pi, 3),
            "radiating_dof": list(body.dofs),
            "water_depth": [8.0],
            "rho": [1025],
            })

    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(
        test_matrix,
        symmetric_body,
        n_threads=1,  # No parallel resolution for cleaner benchmark
        progress_bar=False,
    )

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
