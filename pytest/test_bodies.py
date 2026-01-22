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
    assert np.allclose(body.dofs["1"], np.array([1.0, 0.0, 0.0]))

    body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    assert np.allclose(body.dofs["2"], np.array([0.0, 1.0, 0.0]))

    body.add_rotation_dof(direction=(0.0, 0.0, 1.0), name="3")
    body.add_rotation_dof(rotation_center=(0.5, 0, 0), direction=(0.0, 0.0, 1.0), name="4")


def test_dof_name_inference():
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder())
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_1'])

    body.add_rotation_dof(name="Pitch")
    body.add_rotation_dof(name="yaw")

    body.dofs.clear()
    body.add_all_rigid_body_dofs()


def test_rigid_body_dofs():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(10.0, 0.0, 0.0)))
    assert "Heave" in body.dofs
    assert np.all(body.dofs["Pitch"][:, 2] > 9.0)


def test_rigid_body_dofs_no_rotation_center_but_a_center_of_mass():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(), center_of_mass=(-10.0, 0.0, 0.0))
    assert np.all(body.dofs["Pitch"][:, 2] < -9.0)


def test_rigid_body_dofs_both_a_rotation_center_and_a_center_of_mass():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(10.0, 0.0, 0.0)),
                            center_of_mass=(-10.0, 0.0, 0.0))
    assert np.all(body.dofs["Pitch"][:, 2] > 9.0)


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
    assert body.dofs["Heave"].shape[0] == body.mesh.nb_faces == 1


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
        assert np.allclose(clipped_sphere.dofs['test_dof'], other_clipped_sphere.dofs['test_dof'])
    else:
        assert len(clipped_sphere.dofs['test_dof']) == 0


def test_complicated_clipping_of_dofs():
    # 1 face becomes 2 faces after clipping
    mesh = cpt.Mesh(vertices=[[0.0, 0.0, 0.5], [-0.5, 0.0, -0.5], [0.0, 0.0, -1.5], [0.0, 0.5, -0.5]], faces=[[0, 1, 2, 3]])
    body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())
    clipped_body = body.immersed_part()
    assert len(clipped_body.dofs["Heave"]) == clipped_body.mesh.nb_faces


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
    assert len(clipped_body.dofs["Heave"]) == clipped_body.mesh.nb_faces


def test_cropping_body_with_manual_dof():
    # https://github.com/capytaine/capytaine/issues/204
    sphere = cpt.FloatingBody(cpt.mesh_sphere())
    sphere.dofs["Surge"] = [(1, 0, 0) for face in sphere.mesh.faces]
    sphere = sphere.immersed_part()


def test_immersed_part():
    full_sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(), name="ball")
    immersed_sphere = full_sphere.immersed_part()
    assert immersed_sphere is not full_sphere
    assert immersed_sphere.mesh.nb_faces == full_sphere.mesh.immersed_part().nb_faces
    assert immersed_sphere.mesh.vertices[:, 2].max() <= 0.0


def test_assemble_regular_array():
    body = cpt.FloatingBody(mesh=cpt.mesh_sphere())
    body.add_all_rigid_body_dofs()
    array = body.assemble_regular_array(distance=2.0, nb_bodies=(2, 3))
    assert array.mesh.nb_faces == 6*body.mesh.nb_faces

    # Check name and order of the dofs
    assert list(array.dofs.keys())[0:3] == ["0_0__Surge", "0_0__Sway", "0_0__Heave"]
    assert "2_1__Heave" not in array.dofs.keys()

    # Check that the dofs corresponds to the right panels
    faces_1_0 = np.where(array.dofs["1_0__Heave"] != 0.0)[0]
    fc_1_0 = array.mesh.merged().faces_centers[faces_1_0, :]
    assert np.all(1.0 <= fc_1_0[:, 0]) and np.all(fc_1_0[:, 0] <= 3.0)  #   1 < x < 3
    assert np.all(-1.0 <= fc_1_0[:, 1]) and np.all(fc_1_0[:, 1] <= 1.0) #  -1 < y < 1


r = 1.0
locations = np.array([[1,0],[-1,0],[0,1]])*r*2
n_bodies = locations.shape[0]

@pytest.fixture
def fb_array():
    sphere = cpt.FloatingBody(cpt.mesh_sphere(radius=r, center=(0, 0, 0), resolution=(10, 4)))
    sphere.add_rotation_dof(rotation_center=(0, 0, 0), direction=(0, 1, 0))
    sphere = sphere.immersed_part()
    return sphere.assemble_arbitrary_array(locations)


def test_consistent_dofs_to_faces(fb_array):
    num_active_faces = []
    for fb_dof in fb_array.dofs.items():
        num_active_faces.append(np.count_nonzero(np.count_nonzero(fb_dof[1],axis=1)))

    ma = np.array(num_active_faces)

    tot_faces = fb_array.mesh.nb_faces
    exp_active_faces_per_dof = int(tot_faces/n_bodies)

    assert np.all(ma==exp_active_faces_per_dof)


def test_solve_hydrodynamics(fb_array):
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(coords={
          'rho': 1e3,
          'water_depth': [np.inf],
          'omega': np.pi * 2 / 1,
          'wave_direction': 0,
          'radiating_dof': list(fb_array.dofs.keys()),
          })
    data = solver.fill_dataset(test_matrix, [fb_array],
                                 mesh=True,
                                 wavelength=True,
                                 wavenumber=True)
    assert data.influenced_dof.size == n_bodies
    assert data.radiating_dof.size == n_bodies
    assert data.added_mass.notnull().all()
    assert data.radiation_damping.notnull().all()
    assert data.diffraction_force.notnull().all()
    assert data.Froude_Krylov_force.notnull().all()


# Outdated by v3
def test_clip_component_of_multibody():
    # https://github.com/capytaine/capytaine/issues/660
    body_1 = cpt.FloatingBody(
        mesh=cpt.mesh_sphere(center=(0.0, 0.0, 0.0)),
        dofs=cpt.rigid_body_dofs(rotation_center=(0.0, 0.0, 0.0)),
        name="body_1"
    )
    body_2 = cpt.FloatingBody(
        mesh=cpt.mesh_sphere(center=(5.0, 0.0, 0.0)),
        dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 0.0, 0.0)),
        name="body_2"
    )
    both = body_1 + body_2
    body_2 = body_2.immersed_part()
    assert both.dofs["body_1__Heave"].shape[0] == both.mesh.nb_faces
