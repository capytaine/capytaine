import pytest
import numpy as np
from numpy.linalg import norm
import capytaine as cpt


sphere = cpt.mesh_sphere(radius=1)


def test_clipper():
    """Test clipping of mesh."""
    mesh = cpt.mesh_sphere(radius=5.0, resolution=(10, 5))
    aabb = mesh.axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, water_depth=np.inf)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, water_depth=1.0)
    assert np.allclose(mesh.axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed

    # With CollectionOfMeshes (AxialSymmetry)
    mesh = cpt.mesh_sphere(radius=5.0, resolution=(10, 5), axial_symmetry=True)
    aabb = mesh.merged().axis_aligned_bbox

    mesh.keep_immersed_part(free_surface=0.0, water_depth=np.inf)
    assert np.allclose(mesh.merged().axis_aligned_bbox, aabb[:5] + (0,))  # the last item of the tuple has changed

    mesh.keep_immersed_part(free_surface=0.0, water_depth=1.0)
    assert np.allclose(mesh.merged().axis_aligned_bbox, aabb[:4] + (-1, 0,))  # the last item of the tuple has changed

    # Check boundaries after clipping
    mesh = cpt.mesh_rectangle(size=(5,5), normal=(1,0,0))
    assert max([i[2] for i in mesh.immersed_part(free_surface=-1).vertices])<=-1
    assert max([i[2] for i in mesh.immersed_part(free_surface= 1).vertices])<= 1
    assert min([i[2] for i in mesh.immersed_part(free_surface=100, sea_bottom=-1).vertices])>=-1
    assert min([i[2] for i in mesh.immersed_part(free_surface=100, sea_bottom= 1).vertices])>= 1

    mesh = cpt.mesh_rectangle(size=(4,4), resolution=(1,1), normal=(1,0,0))
    tmp = list(mesh.clip(cpt.Plane(normal=(0,0.1,1),point=(0,0,-1)),inplace=False).vertices)
    tmp.sort(key=lambda x: x[2])
    tmp.sort(key=lambda x: x[1])
    assert np.allclose([i[2] for i in tmp], [-2, -0.8, -2, -1.2])


@pytest.mark.parametrize("size", [5, 6])
def test_clipper_indices(size):
    """Test clipped_mesh_faces_ids."""
    mesh = cpt.mesh_rectangle(size=(size, size), resolution=(size, size), center=(0, 0, 0))
    clipped_mesh = mesh.clipped(plane=cpt.Plane(point=(0, 0, 0), normal=(0, 0, 1)))
    faces_ids = clipped_mesh._clipping_data['faces_ids']

    assert clipped_mesh.nb_faces == len(faces_ids)
    assert all(norm(clipped_mesh.faces_centers[i] - mesh.faces_centers[face_id]) < 0.3
               for i, face_id in enumerate(faces_ids))


def test_clipper_corner_cases():
    mesh = sphere.translated_z(10.0)

    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, 1))
    clipped_mesh = mesh.clip(plane, inplace=False)
    assert clipped_mesh == cpt.Mesh(None, None)  # Empty mesh

    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, -1))
    clipped_mesh = mesh.clip(plane, inplace=False)
    assert clipped_mesh == mesh  # Unchanged mesh

    # Two distinct bodies
    two_spheres = cpt.Mesh.join_meshes(sphere.translated_z(10.0), sphere.translated_z(-10.0))
    plane = cpt.Plane(point=(0, 0, 0), normal=(0, 0, -1))
    one_sphere_remaining = two_spheres.clip(plane, inplace=False)
    assert one_sphere_remaining == sphere.translated_z(10.0)


def test_clipper_tolerance():
    mesh = cpt.mesh_vertical_cylinder(length=10.001, center=(0, 0, -5))
    mesh = mesh.immersed_part()
    np.testing.assert_allclose(mesh.vertices[:, 2].max(), 0.0, atol=1e-12)


def test_degenerate_faces():
    vertices = np.array([
        [-8.00000000e+00,  1.65358984e+00, -4.99999996e-02],
        [-8.00000000e+00,  1.65358984e+00,  5.00000003e-02],
        [-8.00000000e+00,  1.74019238e+00, -9.99999998e-02],
        [-8.00000000e+00,  1.74019238e+00, -1.78037182e-10],
        [-8.00000000e+00,  1.74019238e+00,  1.00000000e-01],
        [-8.00000000e+00,  1.82679492e+00, -5.00000002e-02],
        [-8.00000000e+00,  1.82679492e+00,  4.99999997e-02]
        ])
    faces = np.array([
        [5, 3, 6, 5],
        [5, 2, 3, 5],
        [3, 4, 6, 3],
        [3, 4, 6, 3],
        [2, 0, 3, 2],
        [0, 1, 3, 0],
        [1, 4, 3, 1]
        ])
    mesh = cpt.Mesh(vertices, faces)
    clipped_mesh = mesh.immersed_part()
    assert clipped_mesh.nb_faces == 4
