import logging
import numpy as np
from numpy import pi
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

def make_cylinder(resolution):
    """Make cylinder with a mesh of a given resolution in panels/meter."""
    radius = 1.0
    length = 5.0
    body = cpt.HorizontalCylinder(
        length=length, radius=radius,
        center=(0, 0, -1.5*radius),
        nr=int(resolution*radius),
        nx=int(length*resolution),
        ntheta=int(2*pi*length*resolution),
    )
    # Store the number of panels in the name of the body
    body.name = f"cylinder_{body.mesh.nb_faces:04d}"
    body.add_translation_dof(name="Heave")
    return body

test_matrix = xr.Dataset(coords={
    'omega': [1.0],
    'radiating_dof': ['Heave'],
})

bodies = [make_cylinder(n) for n in np.linspace(1, 5, 10)]
ds1 = cpt.BEMSolver(green_function=cpt.XieDelhommeau()).fill_dataset(test_matrix, bodies)

def read_nb_faces_in_mesh_name(ds):
    """Read the name of the body to guess the resolution of the mesh."""
    ds.coords['nb_faces'] = xr.DataArray([int(name[9:]) for name in ds['body_name'].values], coords=[ds['body_name']])
    ds = ds.swap_dims({'body_name': 'nb_faces'})
    return ds
ds1 = read_nb_faces_in_mesh_name(ds1)

ds1['added_mass'].plot(x='nb_faces')

plt.show()
