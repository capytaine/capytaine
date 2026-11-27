import numpy as np
import xarray as xr

import capytaine as cpt

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
nf_mdf = cpt.compute_near_field_mean_drift_force(dataset, body, rao=rao, solver=solver)
ff_mdf = cpt.compute_far_field_mean_drift_force(dataset, rao=rao)

with np.printoptions(precision=2, suppress=True):
    print('Surge')
    print("Far field: ", ff_mdf['mean_drift_force_far_field'].sel(influenced_dof='Surge').values)
    print("Near field: ", nf_mdf['mean_drift_force_near_field'].sel(influenced_dof='Surge').values)

    print()
    print('Yaw')
    print("Far field: ", ff_mdf['mean_drift_force_far_field'].sel(influenced_dof='Yaw').values)
    print("Near field: ", nf_mdf['mean_drift_force_near_field'].sel(influenced_dof='Yaw').values)
