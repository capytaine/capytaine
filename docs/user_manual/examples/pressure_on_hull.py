import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

mesh = cpt.mesh_sphere(faces_max_radius=0.1).immersed_part()
body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)))

pb = cpt.RadiationProblem(wavelength=1.0, body=body, radiating_dof="Surge", water_depth=10.0)
solver = cpt.BEMSolver()
res = solver.solve(pb, keep_details=True)

body.show_matplotlib(
        color_field=np.real(res.pressure),
        cmap=plt.get_cmap("viridis"),  # Colormap
        )
