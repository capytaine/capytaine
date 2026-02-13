import numpy as np
from numpy import pi
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

cpt.set_logging("INFO")


def make_cylinder(resolution):
    """Make cylinder with a mesh of a given resolution in panels/meter."""
    radius = 1.0
    length = 5.0
    mesh = cpt.mesh_horizontal_cylinder(
        length=length,
        radius=radius,
        center=(0, 0, 0),
        resolution=(
            int(resolution * radius),
            int(2 * pi * radius * resolution),
            int(length * resolution),
        ),
    ).immersed_part()
    body = cpt.FloatingBody(
        mesh,
        lid_mesh=mesh.generate_lid(),
        dofs=cpt.rigid_body_dofs(only=["Heave"]),
        name=f"cylinder_{mesh.nb_faces:04d}"
    )
    # Store the number of panels in the name of the body
    return body


test_matrix = xr.Dataset(
    coords={
        "wavelength": [2.0],
        "radiating_dof": ["Heave"],
    }
)

bodies = [make_cylinder(n) for n in np.linspace(1.0, 9.0, 8)]
ds1 = cpt.BEMSolver().fill_dataset(test_matrix, bodies)


def read_nb_faces_in_mesh_name(ds):
    """Read the name of the body to guess the resolution of the mesh."""
    ds.coords["nb_faces"] = xr.DataArray(
        [int(name[9:]) for name in ds["body_name"].values], coords=[ds["body_name"]]
    )
    ds = ds.swap_dims({"body_name": "nb_faces"})
    return ds


ds1 = read_nb_faces_in_mesh_name(ds1)

ds1["added_mass"].plot(x="nb_faces")

plt.show()
