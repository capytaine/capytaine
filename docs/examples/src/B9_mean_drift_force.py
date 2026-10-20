import numpy as np
import xarray as xr

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force, near_field_mean_drift_force

mesh = cpt.mesh_parallelepiped(resolution=(30, 30, 30)).immersed_part()
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(0,0,0))

wave_direction = [0, np.pi/4]
test_matrix = xr.Dataset(coords={
    'wavenumber': [1.0, 2.0],
    'wave_direction': wave_direction,
    'theta': np.linspace(-np.pi/8, 2*np.pi, 100),  # For far-field mean-drift-force.
    'radiating_dof': list(body.dofs.keys())
})

solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, body, keep_details=True)
rao = cpt.post_pro.rao(dataset)
nf_mdf = near_field_mean_drift_force(rao, dataset, body, solver=solver)
ff_mdf = far_field_mean_drift_force(rao, dataset)

with np.printoptions(precision=2, suppress=True):
    print('Surge')
    print("Far field: ", ff_mdf['drift_force_surge'].values)
    print("Near field: ", nf_mdf['near_field_mean_drift_force'].sel(influenced_dof='Surge').values)

    print()
    print('Yaw')
    print("Far field: ", ff_mdf['drift_force_yaw'].values)
    print("Near field: ", nf_mdf['near_field_mean_drift_force'].sel(influenced_dof='Yaw').values)
