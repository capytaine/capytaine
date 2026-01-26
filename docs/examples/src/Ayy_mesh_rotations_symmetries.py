import numpy as np
import capytaine as cpt
import xarray as xr
from capytaine.new_meshes.symmetric_meshes import RotationSymmetricMesh

cpt.set_logging('WARNING')

meridian_points = np.array([(np.sqrt(1-z**2), 0.0, z) for z in np.linspace(-1.0, 1.0, 50)])
sphere_mesh = RotationSymmetricMesh.from_profile_points(meridian_points, n=50).immersed_part()


symmetric_body = cpt.FloatingBody(
        mesh=sphere_mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        name='symmetric'
        )

full_body = cpt.FloatingBody(
        mesh=sphere_mesh.merged(),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        name='full'
        )

for body in [full_body, symmetric_body]:
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
        body,
        n_threads=1,  # No parallel resolution for cleaner benchmark
        progress_bar=False,
    )

    print("--->", body.name)
    print("Added mass", dataset.added_mass.sel(radiating_dof='Surge', influenced_dof='Surge').values)
    print("Diffraction", dataset.diffraction_force.real.sel(wave_direction=0.0, influenced_dof='Surge').values)
    print(solver.timer_summary())
    print()
