import capytaine as cpt
import numpy as np

def stiffness(mesh):
    floating_body = cpt.FloatingBody(
        mesh,
        dofs=cpt.rigid_body_dofs(only=['Pitch']),
        center_of_mass=(0, 0, 0)
    )
    return floating_body.compute_hydrostatic_stiffness().values

resolutions = [0.6, 0.3, 0.2, 0.1, 0.05, 0.025, 0.01]
# meshes = [cpt.mesh_horizontal_cylinder(length=1, radius=0.5, faces_max_radius=r).immersed_part() for r in resolutions]
meshes = [cpt.mesh_parallelepiped(size=(1, 1, 1), faces_max_radius=r) for r in resolutions]
nb_faces = [mesh.nb_faces for mesh in meshes]
stiffnesses = [stiffness(mesh) for mesh in meshes]
h55_first_order = np.array([h[0, 0] for h in stiffnesses])

stiffnesses_with_quad = [stiffness(mesh.with_quadrature('Gauss-Legendre 2')) for mesh in meshes]
h55_quad = np.array([h[0, 0] for h in stiffnesses_with_quad])

import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].plot(nb_faces, np.abs(h55_first_order), label='first order quadrature')
axs[0].plot(nb_faces, np.abs(h55_quad), label='Gauss-Legendre 2 quadrature')
axs[1].plot(nb_faces[:-1], np.abs(h55_first_order[:-1] - h55_quad[-1]), linestyle='--', label='erreur, first order quadrature')
axs[1].plot(nb_faces[:-1], np.abs(h55_quad[:-1] - h55_quad[-1]), linestyle='--', label='erreur, Gauss-Legendre 2 quadrature')
for ax in axs:
    ax.set_xlabel('resolution')
    ax.set_xscale('log')
    ax.set_yscale('log')
fig.suptitle('Hydrostatic stiffness h55 vs nb faces')
plt.legend()
plt.show()
