import capytaine as cpt
import numpy as np
import xarray as xr
import timeit
import os
import matplotlib.pyplot as plt

cpu_count = os.cpu_count()
print(cpu_count)

sphere = cpt.FloatingBody(
        cpt.mesh_sphere(radius=1, center=(0, 0, -2.0),
                        resolution=(20, 20)),
        name="sphere")
sphere.add_translation_dof(name="Heave")

test_matrix = xr.Dataset(coords={
    'period': np.linspace(1, 3, 200),
    'wave_direction': [0],
    'radiating_dof': list(sphere.dofs),
    'water_depth': [np.inf]})
solver = cpt.BEMSolver()

n_workers = [1] + [int(n*2) for n in np.arange(1,np.floor(cpu_count/2)+1)]
times = np.array([timeit.timeit(
    lambda : solver.fill_dataset(test_matrix, sphere, n_jobs=n_jobs), 
    number=1, # to be more rigorous, this should be greater than 1
    ) 
                  for n_jobs in n_workers])

xx = np.linspace(1, cpu_count)

fig, ax = plt.subplots()
ax.plot(n_workers, times[0]/times, marker='.', 
        label='Capytaine on your machine')
ax.plot(xx,xx,'--', label='Linear')
ax.set_xlabel('Cores [ ]')
ax.set_ylabel('Speedup [ ]')
ax.legend()
plt.show()
