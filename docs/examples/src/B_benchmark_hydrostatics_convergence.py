import capytaine as cpt
import numpy as np

def h55_stiffness(mesh):
    floating_body = cpt.FloatingBody(
        mesh,
        dofs=cpt.rigid_body_dofs(only=['Pitch'], rotation_center=(0.25, 0.25, -0.25)),
        center_of_mass=(0.25, 0.25, -0.25)
    )
    return floating_body.compute_hydrostatic_stiffness().values[0, 0]

def volume(mesh):
    return mesh.volume

resolutions = [0.2, 0.1, 0.05, 0.025, 0.01]
# meshes = [cpt.mesh_horizontal_cylinder(length=1, radius=0.5, faces_max_radius=r).immersed_part() for r in resolutions]
meshes = [cpt.mesh_parallelepiped(size=(1, 1, 1), faces_max_radius=r).immersed_part() for r in resolutions]

meshes_quad = [mesh.with_quadrature('Gauss-Legendre 2') for mesh in meshes]
nb_faces = [mesh.nb_faces for mesh in meshes]
nb_faces_range = np.array(nb_faces[:-1])

stiffnesses = np.array([h55_stiffness(mesh) for mesh in meshes])
stiffnesses_with_quad = np.array([h55_stiffness(mesh) for mesh in meshes_quad])

volumes = np.array([volume(mesh) for mesh in meshes])
volumes_quad = np.array([volume(mesh) for mesh in meshes_quad])

import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 2, sharex=True, layout="constrained")

# Volume: value
axs[0, 0].plot(nb_faces, volumes, label='first order quadrature')
axs[0, 0].plot(nb_faces, volumes_quad, label='Gauss-Legendre 2 quadrature')
axs[0, 0].set_ylabel('volume')
axs[0, 0].set_title('Volume')

# Volume: error
axs[1, 0].plot(nb_faces_range, np.abs(volumes[:-1] - volumes_quad[-1]), linestyle='--', label='erreur, first order quadrature')
axs[1, 0].plot(nb_faces_range, np.abs(volumes_quad[:-1] - volumes_quad[-1]), linestyle='--', label='erreur, Gauss-Legendre 2 quadrature')
ref = np.abs(volumes_quad[:-1] - volumes_quad[-1])[0] * (nb_faces_range[0] / nb_faces_range) ** 1
axs[1, 0].plot(nb_faces_range, ref, 'k:', alpha=0.5, linewidth=1, label='O(h²) = O(N⁻¹)')
axs[1, 0].set_ylabel('error')

# Stiffness h55: value
axs[0, 1].plot(nb_faces, np.abs(stiffnesses), label='first order quadrature')
axs[0, 1].plot(nb_faces, np.abs(stiffnesses_with_quad), label='Gauss-Legendre 2 quadrature')
axs[0, 1].set_ylabel('|h55|')
axs[0, 1].set_title('Stiffness h55')

# Stiffness h55: error
axs[1, 1].plot(nb_faces_range, np.abs(stiffnesses[:-1] - stiffnesses_with_quad[-1]), linestyle='--', label='erreur, first order quadrature')
axs[1, 1].plot(nb_faces_range, np.abs(stiffnesses_with_quad[:-1] - stiffnesses_with_quad[-1]), linestyle='--', label='erreur, Gauss-Legendre 2 quadrature')
ref = np.abs(stiffnesses_with_quad[:-1] - stiffnesses_with_quad[-1])[0] * (nb_faces_range[0] / nb_faces_range) ** 1
axs[1, 1].plot(nb_faces_range, ref, 'k:', alpha=0.5, linewidth=1, label='O(h²) = O(N⁻¹)')
axs[1, 1].set_ylabel('error')


for ax in axs.flat:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize='small')
for ax in axs[1, :]:
    ax.set_xlabel('nb_faces')
fig.suptitle('Hydrostatic quantities vs nb faces')
plt.show()
