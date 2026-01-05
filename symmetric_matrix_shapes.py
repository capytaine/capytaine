
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

from capytaine.new_meshes import (
        Mesh,
    RotationSymmetricMesh,
    ReflectionSymmetricMesh
)
from capytaine.new_meshes.meshes import to_new_mesh

# reference_mesh = cpt.mesh_parallelepiped(resolution=(8, 12, 9), name="reference_mesh")
# half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
# quarter_mesh = half_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(0, 1, 0)))
# sym_mesh = ReflectionSymmetricMesh(half=to_new_mesh(half_mesh), plane="yOz", name="new_simple_symmetric_mesh")
# # sym_mesh = ReflectionSymmetricMesh(half=ReflectionSymmetricMesh(half=to_new_mesh(quarter_mesh), plane="xOz"), plane="yOz", name="new_nested_symmetric_mesh")

reference_mesh = cpt.mesh_sphere(resolution=(10, 10), name="reference_mesh")
new_reference_mesh = to_new_mesh(reference_mesh)
wedge = new_reference_mesh.extract_wedge(n=5)
half_wedge = wedge.rotated_z(-np.pi/5).clipped(origin=(0, 0, 0), normal=(0.0, 1.0, 0.0))
# sym_mesh = RotationSymmetricMesh(wedge=wedge, n=5)
sym_mesh = RotationSymmetricMesh(ReflectionSymmetricMesh(half_wedge, plane="xOz"), n=5)
# sym_mesh = ReflectionSymmetricMesh(RotationSymmetricMesh(wedge=half_wedge, n=5), plane="xOz")

half_sphere = to_new_mesh(cpt.mesh_sphere(center=(1.0, 0.0, 0.0), radius=0.2, resolution=(4, 4)).immersed_part()).clipped(origin=(1, 0, 0), normal=(0, 1, 0))
sym_mesh = ReflectionSymmetricMesh(RotationSymmetricMesh(half_sphere, n=10), plane="xOz")
# sym_mesh = RotationSymmetricMesh(ReflectionSymmetricMesh(half_sphere, plane="xOz"), n=10)
# sym_mesh = ReflectionSymmetricMesh(RotationSymmetricMesh(half_sphere, n=3).merged(), plane="xOz")

ref_mesh = sym_mesh.merged()
class DistanceKernelFunction:
    exportable_settings = {}
    def evaluate(self, mesh1, mesh2, *args, **kwargs):
        A = np.array([[area*np.linalg.norm(x-xi) for xi in mesh2.faces_centers] for (x, area) in zip(mesh1.faces_centers, mesh1.faces_areas)])
        return A, A
engine = cpt.BasicMatrixEngine(
    green_function=DistanceKernelFunction()
)
params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
S_ref, K_ref = engine.build_matrices(ref_mesh, ref_mesh, **params)
S, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
# # S_ref, K_ref = engine.build_matrices(sym_mesh.half.merged(), sym_mesh.other_half.merged(), **params)
# S_ref, K_ref = engine.build_matrices(sym_mesh.half.merged(), sym_mesh.half.merged().mirrored('xOz'), **params)
# S, K = engine.build_matrices(sym_mesh.half, sym_mesh.other_half, **params)
np.linalg.norm(sym_mesh.faces_centers - sym_mesh.merged().faces_centers)
np.linalg.norm(sym_mesh.mirrored("xOz").faces_centers - sym_mesh.merged().mirrored("xOz").faces_centers)
fig, axs = plt.subplots(3, 1)
ims = axs[0].imshow(np.real(S_ref))
axs[0].set_title('ref')
plt.colorbar(ims, ax=axs[0])
ims = axs[1].imshow(np.real(np.array(S)))
axs[1].set_title('sym')
plt.colorbar(ims, ax=axs[1])
ims = axs[2].imshow(np.abs(S_ref - np.array(S)), cmap='magma', vmin=0.0)
axs[2].set_title('diff')
plt.colorbar(ims, ax=axs[2])
plt.show()
