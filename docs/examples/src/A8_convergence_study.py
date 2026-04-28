import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

cpt.set_logging("INFO")

resolutions = np.linspace(1.0, 8.0, 6)

test_matrix = xr.Dataset(
    coords={
        "wavelength": [2.0],
        "radiating_dof": ["Heave"],
    }
)
solver = cpt.BEMSolver()

datasets = []
for resolution in resolutions:
    radius = 1.0
    length = 5.0

    mesh = cpt.mesh_horizontal_cylinder(
        length=length,
        radius=radius,
        center=(0, 0, 0),
        resolution=(
            int(resolution * radius),
            int(2 * np.pi * radius * resolution),
            int(length * resolution),
        ),
    ).immersed_part()

    body = cpt.FloatingBody(
        mesh,
        lid_mesh=mesh.generate_lid(),
        dofs=cpt.rigid_body_dofs(only=["Heave"]),
        name=f"cylinder_with_resolution_{resolution:.1f}"
    )

    datasets.append(solver.fill_dataset(test_matrix, body, mesh=True))
    # When `mesh=True`, some extra data about the mesh is stored in the dataset, including the number of faces.


nb_faces = [ds.coords["nb_faces"] for ds in datasets]
added_mass = [ds["added_mass"].values.squeeze() for ds in datasets]
plt.plot(nb_faces, added_mass)
plt.xlabel("nb faces")
plt.ylabel("added mass")
plt.show()
