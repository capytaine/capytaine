
from capytaine.geometric_bodies.sphere import Sphere
from capytaine.geometric_bodies.cylinder import HorizontalCylinder
from capytaine.problems import DiffractionProblem, RadiationProblem
from capytaine.results import assemble_dataset, wavenumber_data_array
from capytaine.Nemoh import Nemoh

sphere = Sphere(radius=1.0, center=(0, 0, -2.0),
                ntheta=20, nphi=20,
                name="sphere")
sphere.add_translation_dof(name="Surge")
sphere.add_translation_dof(name="Heave")

cylinder = HorizontalCylinder(length=5.0, radius=1.0,
                              center=(+1.5, 3.0, -3.0),
                              nx=20, nr=3, ntheta=20,
                              name="cylinder")
cylinder.add_translation_dof(name="Surge")
cylinder.add_translation_dof(name="Heave")

both = cylinder + sphere
both.show()

solver = Nemoh()
problems = [RadiationProblem(body=both, radiating_dof=dof, omega=1.0) for dof in both.dofs]
problems += [DiffractionProblem(body=both, angle=0.0, omega=1.0)]
results = solver.solve_all(problems)
data = assemble_dataset(results)

print(data)
