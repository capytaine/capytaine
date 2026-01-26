import pytest

import numpy as np
import xarray as xr

import capytaine as cpt


def test_dof():
    nodes = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]])
    faces = np.array([[0, 1, 2, 3]])
    body = cpt.FloatingBody(cpt.Mesh(nodes, faces), name="one_face")
    assert body.dofs == {}

    body.add_translation_dof(direction=(1.0, 0.0, 0.0), name="1")
    assert np.allclose(
        body.dofs["1"].evaluate_motion(body.mesh),
        np.array([1.0, 0.0, 0.0])
    )

    body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    assert np.allclose(
        body.dofs["2"].evaluate_motion(body.mesh),
        np.array([0.0, 1.0, 0.0])
    )

    body.add_rotation_dof(direction=(0.0, 0.0, 1.0), name="3")
    body.add_rotation_dof(rotation_center=(0.5, 0, 0), direction=(0.0, 0.0, 1.0), name="4")


def test_dof_name_inference():
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder())
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(
            body.dofs[dofname].evaluate_motion(body.mesh),
            body.dofs['Surge_1'].evaluate_motion(body.mesh)
        )

    body.add_rotation_dof(name="Pitch")
    body.add_rotation_dof(name="yaw")

    body.dofs.clear()
    body.add_all_rigid_body_dofs()


def test_rigid_body_dofs():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(10.0, 0.0, 0.0)))
    assert "Heave" in body.dofs
    assert np.all(body.dofs["Pitch"].evaluate_motion(mesh)[:, 2] > 9.0)


def test_rigid_body_dofs_both_a_rotation_center_and_a_center_of_mass():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(10.0, 0.0, 0.0)),
                            center_of_mass=(-10.0, 0.0, 0.0))
    assert np.all(body.dofs["Pitch"].evaluate_motion(mesh)[:, 2] > 9.0)


def test_rigid_body_dofs_neither_a_rotation_center_nor_a_center_of_mass():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    assert np.allclose(body._infer_rotation_center(), (0.0, 0.0, 0.0))


def test_defining_rotation_center_with_ints():
    # Issue #319
    mesh = cpt.mesh_sphere().immersed_part()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -1)))
    body.translated_y(-2.0)


def test_healing_before_initializing_dofs():
    # Issue #367: https://github.com/capytaine/capytaine/issues/367
    # Define a mesh with a normal panel and degenerate panel
    vertices = np.array([(0.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                         (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)])
    faces = np.array([[0, 1, 2, 3], [1, 2, 2, 1]])
    mesh = cpt.Mesh(vertices, faces)
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    assert body.dofs["Heave"].evaluate_motion(mesh).shape[0] == body.mesh.nb_faces == 1


def test_translate_center_of_mass_defined_as_tuple():
    body = cpt.FloatingBody(center_of_mass=(0, 0, 0))
    body = body.translated((0, 0, 1))
    assert body.center_of_mass.shape == (3,)
    np.testing.assert_allclose(body.center_of_mass, np.array([0, 0, 1]))


def test_rotate_center_of_mass_defined_as_tuple():
    body = cpt.FloatingBody(center_of_mass=(1, 0, 0))
    body = body.rotated_z(np.pi)
    assert body.center_of_mass.shape == (3,)
    np.testing.assert_allclose(body.center_of_mass, np.array([-1, 0, 0]), atol=1e-10)


def test_mirror_center_of_mass_defined_as_tuple():
    body = cpt.FloatingBody(center_of_mass=(0, -1, 0))
    body = body.mirrored('xOz')
    assert body.center_of_mass.shape == (3,)
    np.testing.assert_allclose(body.center_of_mass, np.array([0, -1, 0]))


def test_transform_with_lid():
    mesh = cpt.mesh_sphere().immersed_part()
    lid_mesh = mesh.generate_lid()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh)
    body.mesh_including_lid
    body = body.translated_x(10.0)
    assert np.all(body.mesh.vertices[:, 0] > 5.0)
    assert np.all(body.mesh_including_lid.vertices[:, 0] > 5.0)


@pytest.mark.parametrize("z_center", [0, 2, -2])
@pytest.mark.parametrize("symmetry", [True, False])
def test_clipping_of_dofs(z_center, symmetry):
    """Check that clipping a body with a dof is the same as clipping the body ant then adding the dof."""
    full_sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(center=(0, 0, z_center), axial_symmetry=symmetry), name="sphere")

    full_sphere.add_rotation_dof(rotation_center=(1, 0, 0), direction=(1, 0, 0), name="test_dof")
    clipped_sphere = full_sphere.immersed_part(free_surface=0.0, water_depth=np.inf)

    other_clipped_sphere = cpt.FloatingBody(mesh=clipped_sphere.mesh, name="other_sphere")
    other_clipped_sphere.add_rotation_dof(rotation_center=(1, 0, 0), direction=(1, 0, 0), name="test_dof")

    if clipped_sphere.mesh.nb_faces > 0:
        assert np.allclose(
            clipped_sphere.dofs['test_dof'].evaluate_motion(clipped_sphere.mesh),
            other_clipped_sphere.dofs['test_dof'].evaluate_motion(other_clipped_sphere.mesh)
        )
    else:
        assert len(clipped_sphere.dofs['test_dof'].evaluate_motion(clipped_sphere.mesh)) == 0


def test_complicated_clipping_of_dofs():
    # 1 face becomes 2 faces after clipping
    mesh = cpt.Mesh(vertices=[[0.0, 0.0, 0.5], [-0.5, 0.0, -0.5], [0.0, 0.0, -1.5], [0.0, 0.5, -0.5]], faces=[[0, 1, 2, 3]])
    body = cpt.FloatingBody(mesh, dofs={'Shear': np.array([(x, 0, 0) for (x, y, z) in mesh.faces_centers])})
    clipped_body = body.immersed_part()
    assert len(clipped_body.dofs["Shear"]) == clipped_body.mesh.nb_faces


def test_clipping_of_dofs_with_degenerate_faces():
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
        [3, 4, 6, 3],
        [2, 0, 3, 2],
        ])
    mesh = cpt.Mesh(vertices, faces)
    body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())
    clipped_body = body.immersed_part()
    assert len(clipped_body.dofs["Heave"].evaluate_motion(clipped_body.mesh)) == clipped_body.mesh.nb_faces


def test_cropping_body_with_manual_dof():
    # https://github.com/capytaine/capytaine/issues/204
    sphere = cpt.FloatingBody(cpt.mesh_sphere())
    sphere.dofs["Shear"] = [(x, 0, 0) for (x, y, z) in sphere.mesh.faces_centers]
    sphere = sphere.immersed_part()


def test_immersed_part():
    full_sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(), name="ball")
    immersed_sphere = full_sphere.immersed_part()
    assert immersed_sphere is not full_sphere
    assert immersed_sphere.mesh.nb_faces == full_sphere.mesh.immersed_part().nb_faces
    assert immersed_sphere.mesh.vertices[:, 2].max() <= 0.0
