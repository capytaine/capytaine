#! /usr/bin/env python3

import numpy as np
import capytaine as cpt


r = 1.
draft = 0.5
depth = 10.
omega = 1.
rho = 1000

# Define geometry and heave degree of freedom
body = cpt.FloatingBody(
        mesh=cpt.mesh_vertical_cylinder(
            length=2*draft, radius=r, center=(0.,0.,0.),
            resolution=(10, 50, 10)
            ))
body.add_translation_dof(name='Heave')
body = body.immersed_part()

solver = cpt.BEMSolver()

# Define and solve the diffraction and radiation problems
diff_problem = cpt.DiffractionProblem(body=body, water_depth=depth,
                                      omega=omega, wave_direction=0.)
rad_problem = cpt.RadiationProblem(body=body, water_depth=depth,
                                   omega=omega, radiating_dof='Heave')

diff_solution = solver.solve(diff_problem)
rad_solution = solver.solve(rad_problem)


# Read mesh properties
faces_centers = body.mesh.faces_centers
faces_normals = body.mesh.faces_normals
faces_areas = body.mesh.faces_areas

from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force

# Computation from the diffraction solution (Capytaine)
FK = froude_krylov_force(diff_problem)['Heave']
diff_force = diff_solution.forces['Heave']

# Get potentials
phi_inc = airy_waves_potential(faces_centers, diff_problem)
v_inc = airy_waves_velocity(faces_centers, diff_problem)
phi_rad = rad_solution.potential

# Indirect computation from the radiation solution, via the Haskind relation
integrand = - (phi_inc * faces_normals[:,2]
              - 1j/omega * phi_rad * np.diag(v_inc@faces_normals.T))
F_hask = 1j * omega * rho * np.sum(integrand*faces_areas)

# Direct integration of the potentials
diff_force_recomputed = - 1j * omega * rho * np.sum(
                        diff_solution.potential*faces_areas*faces_normals[:,2])
FK_recomputed = - 1j * omega * rho * np.sum(
                        phi_inc*faces_areas*faces_normals[:,2])

print('Result from Capytaine = ', FK + diff_force)
print('Result from recomputed direct integration',
        FK_recomputed + diff_force_recomputed)
print('Result from Haskind relation = ', F_hask)
