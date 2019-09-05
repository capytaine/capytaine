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
    radius = 1.0
    length = 3.0
    body = cpt.HorizontalCylinder(
        length=length, radius=radius,
        center=(0, 0, -1.01*radius),
        nr=int(resolution*radius),
        nx=int(length*resolution),
        ntheta=int(2*pi*length*resolution),
    )
    body.name = f"cylinder_{resolution}"
    body.add_translation_dof(name="Heave")
    return body

test_matrix = xr.Dataset(coords={
    'omega': [1.0],
    'radiating_dof': ['Heave'],
})

bodies = [make_cylinder(n) for n in range(3, 10, 2)]

def extract_resolution_from_body_name(ds):
    ds.coords['resolution'] = xr.DataArray([int(name[9:]) for name in ds['body_name'].values], coords=[ds['body_name']])
    ds = ds.swap_dims({'body_name': 'resolution'})
    return ds

ds1 = cpt.BEMSolver(green_function=cpt.XieDelhommeau()).fill_dataset(test_matrix, bodies)
ds1 = extract_resolution_from_body_name(ds1)

plt.figure(1)
ds1['added_mass'].plot(x='resolution', label='Delhommeau')
plt.legend()
plt.figure(2)
ds1['radiation_damping'].plot(x='resolution', label='Delhommeau')
plt.legend()

plt.show()
