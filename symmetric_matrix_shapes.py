
import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt

from capytaine.new_meshes import (
    RotationSymmetricMesh,
    ReflectionSymmetricMesh
)
from capytaine.new_meshes.meshes import to_new_mesh

class DistanceKernelFunction:
    exportable_settings = {}
    def evaluate(self, mesh1, mesh2, *args, **kwargs):
        A = np.array([[area*np.linalg.norm(x-xi) for xi in mesh2.faces_centers] for (x, area) in zip(mesh1.faces_centers, mesh1.faces_areas)])
        return A, A

def reflection():
    reference_mesh = cpt.mesh_parallelepiped(resolution=(8, 12, 9), name="reference_mesh")
    half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
    return ReflectionSymmetricMesh(half=to_new_mesh(half_mesh), plane="yOz", name="new_simple_symmetric_mesh")

def nested_reflections():
    reference_mesh = cpt.mesh_parallelepiped(resolution=(8, 12, 9), name="reference_mesh")
    half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
    quarter_mesh = half_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(0, 1, 0)))
    return ReflectionSymmetricMesh(half=ReflectionSymmetricMesh(half=to_new_mesh(quarter_mesh), plane="xOz"), plane="yOz", name="new_nested_symmetric_mesh")

def rotation(n=4):
    reference_mesh = cpt.mesh_sphere(resolution=(20, 20), name="reference_mesh")
    new_reference_mesh = to_new_mesh(reference_mesh)
    wedge = new_reference_mesh.extract_wedge(n=n)
    return RotationSymmetricMesh(wedge=wedge, n=n)

def rotation_of_reflection(n=4):
    half_sphere = to_new_mesh(cpt.mesh_sphere(center=(1.0, 0.0, 0.0), radius=0.2, resolution=(4, 4)).immersed_part()).clipped(origin=(1, 0, 0), normal=(0, 1, 0))
    return RotationSymmetricMesh(ReflectionSymmetricMesh(half_sphere, plane="xOz"), n=n)

def reflection_of_rotation(n=4):
    half_sphere = to_new_mesh(cpt.mesh_sphere(center=(1.0, 0.0, 0.0), radius=0.2, resolution=(4, 4)).immersed_part()).clipped(origin=(1, 0, 0), normal=(0, 1, 0))
    return ReflectionSymmetricMesh(RotationSymmetricMesh(wedge=half_sphere, n=n), plane="xOz")

def test(sym_mesh):
    ref_mesh = sym_mesh.merged()
    engine = cpt.BasicMatrixEngine(
        green_function=DistanceKernelFunction()
    )
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
    S_ref, K_ref = engine.build_matrices(ref_mesh, ref_mesh, **params)
    S, K = engine.build_matrices(sym_mesh, sym_mesh, **params)

    fig, axs = plt.subplots(2, 2)
    axs = axs.ravel()
    ims = axs[0].imshow(np.real(S_ref))
    axs[0].set_title('ref')
    plt.colorbar(ims, ax=axs[0])
    ims = axs[1].imshow(np.real(np.array(S)))
    axs[1].set_title('sym')
    plt.colorbar(ims, ax=axs[1])
    ims = axs[2].imshow(np.abs(S_ref - np.array(S)), cmap='magma', vmin=0.0)
    axs[2].set_title('diff')
    plt.colorbar(ims, ax=axs[2])

# test(reflection())
# test(nested_reflections())
# test(rotation(n=3))
test(rotation_of_reflection(n=3))
# test(reflection_of_rotation(n=3))

plt.show()
