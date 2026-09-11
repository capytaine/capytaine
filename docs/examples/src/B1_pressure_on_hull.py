import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

mesh = cpt.mesh_sphere(radius=1.0, faces_max_radius=0.1).immersed_part()
body = cpt.FloatingBody(
    mesh=mesh,
    lid_mesh=mesh.generate_lid(),
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0))
)

pb = cpt.RadiationProblem(wavelength=1.0, body=body, radiating_dof="Surge", water_depth=np.inf)
solver = cpt.BEMSolver()
res = solver.solve(pb, keep_details=True)

body.mesh.show(
        backend="matplotlib",
        color_field=np.real(res.pressure_on_hull),
        cmap=plt.get_cmap("viridis"),  # Colormap
        )

# # Plotting also the (not physically meaningful) pressure on the lid
# body.mesh_including_lid.show(
#         backend="matplotlib",
#         color_field=np.real(res.pressure),
#         cmap=plt.get_cmap("viridis"),  # Colormap
#         )
