#! /usr/bin/env python3

import numpy as np
import capytaine as cpt
from time import time

radius = 2.
spacing = 5.
nb_cols = 14
nb_rows = 14

omega = 1.

x, y = np.meshgrid(np.arange(nb_cols)*spacing, np.arange(nb_rows)*spacing)
x = x.flatten()
y = y.flatten()

bodies = [cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=radius,
                                                center=(x[ii], y[ii], 0),
                                                resolution=(10,10)))
            for ii in range(len(x))]

for body in bodies:
    body.keep_immersed_part()
    body.add_translation_dof(name='Heave')

# Dense problem setup and solution
print('Solving dense problem')
all_bodies = cpt.FloatingBody.join_bodies(*bodies, name="joined bodies")
problems = [cpt.RadiationProblem(body=all_bodies, radiating_dof=dof, omega=omega)
        for dof in all_bodies.dofs]
solver = cpt.BEMSolver()
dresults = solver.solve_all(problems)
ddata = cpt.assemble_dataset(dresults)


# Hierarchical problem setup
all_bodies = cpt.FloatingBody.cluster_bodies(*bodies, name="clustered bodies")
problems = [cpt.RadiationProblem(body=all_bodies, radiating_dof=dof, omega=omega)
            for dof in all_bodies.dofs]


cpt.set_logging(level='INFO')  # prints the number of iterations
# Hierarchical solve, with preconditioner
print('Solving hierarchical problem with preconditioner')
t0 = time()
sparse_engine = cpt.HierarchicalPrecondMatrixEngine(ACA_distance=2.0, ACA_tol=1e-2)
solver = cpt.BEMSolver(engine=sparse_engine)
presults = solver.solve_all(problems)
tP = time() - t0

data = cpt.assemble_dataset(presults)

# Hierarchical solve, no preconditioner
print('Solving hierarchical problem without preconditioner')
t0 = time()
sparse_engine = cpt.HierarchicalToeplitzMatrixEngine(ACA_distance=2.0, ACA_tol=1e-2)
solver = cpt.BEMSolver(engine=sparse_engine)
hresults = solver.solve_all(problems)
tNP = time() - t0

# The clustering process changes the order of the bodies with respect to
# the one in list 'bodies'. The correspondence with keys is maintained correct.
# However, if one wants to go back to the initial ordering, the following
# lines can be used.

ordered_dofs_names = []
for body in bodies:
    ordered_dofs_names += [body.name + '__' + dofkey for dofkey in list(body.dofs.keys())  ]
reordered_dofs = {'radiating_dof': ordered_dofs_names,
                  'influenced_dof': ordered_dofs_names}

reordered_data = data.reindex(reordered_dofs)

print('error = ', np.linalg.norm(presults[0].sources - hresults[0].sources))
print('time with preconditioner = ', tP)
print('time without preconditioner = ', tNP)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3)
ax[0].matshow(ddata['added_mass'].data[0,:,:])
ax[1].matshow(reordered_data['added_mass'].data[0,:,:])
ax[2].matshow(data['added_mass'].data[0,:,:])
ax[0].set_title('Dense')
ax[1].set_title('HP, reordered')
ax[2].set_title('HP, before reordering')
fig.suptitle('Added mass matrix')

plt.show()
