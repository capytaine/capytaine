#!/usr/bin/env python
# coding: utf-8

from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.cylinder import HorizontalCylinder
from capytaine.problems import DiffractionProblem, RadiationProblem
from capytaine.Nemoh import Nemoh
from capytaine.results import assemble_dataset

# Define the first body
sphere = Sphere(radius=1.0, center=(0, 0, -2.0),
                ntheta=20, nphi=20,
                name="sphere_1")
sphere.add_translation_dof(name="Surge")
sphere.add_translation_dof(name="Heave")

# Define the second body
other_sphere = Sphere(radius=0.5, center=(-2, -3, -1),
                      ntheta=20, nphi=20,
                      name="sphere_2")
other_sphere.add_translation_dof(name="Surge")
other_sphere.add_translation_dof(name="Heave")

# Define the third body
cylinder = HorizontalCylinder(length=5.0, radius=1.0,
                              center=(1.5, 3.0, -3.0),
                              nx=20, nr=3, ntheta=20,
                              name="cylinder")
cylinder.add_translation_dof(name="Surge")
cylinder.add_translation_dof(name="Heave")

# Combine the three individual bodies into a single body.
all_bodies = cylinder + sphere + other_sphere
#
# Note that the + operation is associative, that is "(cylinder + sphere) + other_sphere"
# is equivalent to "cylinder + (sphere + other_sphere)".
# It is also commutative, up to some internal details which are usually not relevant.

print("Merged body name:", all_bodies.name)
print("Merged body dofs:", list(all_bodies.dofs.keys()))

# The merged body can be used to define the problems in the usual way
problems = [RadiationProblem(body=all_bodies, radiating_dof=dof, omega=1.0) for dof in all_bodies.dofs]
problems += [DiffractionProblem(body=all_bodies, angle=0.0, omega=1.0)]

# Solves the problem
solver = Nemoh()
results = solver.solve_all(problems)
data = assemble_dataset(results)

print(data)
