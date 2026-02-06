import capytaine as cpt
import numpy as np
mesh = cpt.mesh_sphere(radius=1.0)
self = mesh

(mesh.faces_normals[:, None, :] * self.quadrature_points[0] * self.quadrature_points[1])
print(mesh.wet_surface_area)
print(mesh.volumes)
print(mesh.surface_integral(1))
print(mesh.waterplane_integral(1))
print(mesh.center_of_buoyancy)
print(mesh.with_quadrature('Gauss-Legendre 2').center_of_buoyancy)

mesh.with_quadrature('Gauss-Legendre 2').immersed_part()

def stiffness(mesh):
    floating_body = cpt.FloatingBody(
        mesh,
        dofs=cpt.rigid_body_dofs(),
        center_of_mass=(0, 0, 0)
    )
    return floating_body.compute_hydrostatic_stiffness().values

resolutions = [0.3, 0.2, 0.1, 0.05, 0.03]
meshes = [cpt.mesh_horizontal_cylinder(length=1, radius=1, faces_max_radius=r) for r in resolutions]
nb_faces = [mesh.nb_faces for mesh in meshes]
stiffnesses = [stiffness(mesh) for mesh in meshes]
h55_first_order = np.array([h[4, 4] for h in stiffnesses])

stiffnesses_with_quad = [stiffness(mesh.with_quadrature('Gauss-Legendre 2')) for mesh in meshes]
h55_quad = np.array([h[4, 4] for h in stiffnesses_with_quad])

import matplotlib.pyplot as plt
plt.plot(nb_faces, np.abs(h55_first_order), label='default quadrature')
plt.plot(nb_faces, np.abs(h55_quad), label='Gauss-Legendre 2 quadrature')
plt.plot(nb_faces[:-1], np.abs(h55_first_order[:-1] - h55_first_order[-1]), linestyle='--', label='erreur, default quadrature')
plt.plot(nb_faces[:-1], np.abs(h55_quad[:-1] - h55_quad[-1]), linestyle='--', label='erreur, Gauss-Legendre 2 quadrature')
plt.xlabel('resolution')
plt.title('Hydrostatic stiffness h55 vs resolution')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
