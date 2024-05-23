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

    body.add_rotation_dof(cpt.Axis(vector=(0.0, 0.0, 1.0)), name="3")
    body.add_rotation_dof(cpt.Axis(point=(0.5, 0, 0), vector=(0.0, 0.0, 1.0)), name="4")


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


def test_bodies():
    body = cpt.FloatingBody(mesh=cpt.mesh_sphere(), name="sphere", center_of_mass=(0, 0, 0))
    repr(body)
    body.add_translation_dof(name="Surge")
    body.add_translation_dof(name="Heave")

    # Extract faces
    body.extract_faces(np.where(body.mesh.faces_centers[:, 2] < 0)[0])

    # Clipping
    body.keep_immersed_part(inplace=False)

    # Mirror of the dofs
    mirrored = body.mirrored(cpt.Plane(point=(1, 0, 0), normal=(1, 0, 0)))
    assert np.allclose(mirrored.center_of_mass, np.array([2, 0, 0]))
    assert np.allclose(body.dofs['Surge'], -mirrored.dofs['Surge'])

    # Rotation of the dofs
    sideways = body.rotated(cpt.Axis(point=(0, 0, 0), vector=(0, 1, 0)), np.pi/2)
    assert np.allclose(sideways.dofs['Heave'][0], np.array([1, 0, 0]))

    upside_down = body.rotated(cpt.Axis(point=(0, 0, 0), vector=(0, 1, 0)), np.pi)
    assert np.allclose(body.dofs['Heave'], -upside_down.dofs['Heave'])

    # Copy of the body
    copy_of_body = body.copy(name="copy_of_sphere")
    copy_of_body.translate_x(10.0)
    copy_of_body.add_translation_dof(name="Heave")

    # Join bodies
    both = body.join_bodies(copy_of_body)
    assert set(both.dofs) == {'sphere__Surge', 'copy_of_sphere__Surge', 'sphere__Heave', 'copy_of_sphere__Heave'}


def test_mirror_rotation_center_defined_as_tuple():
    body = cpt.FloatingBody()
    body.rotation_center = (0, 1, 0)
    body.mirror(cpt.xOz_Plane)
    assert body.rotation_center.shape == (3,)
    np.testing.assert_allclose(body.rotation_center, np.array([0, -1, 0]))


def test_translate_rotation_center_defined_as_tuple():
    body = cpt.FloatingBody()
    body.rotation_center = (0, 0, 0)
    body.translate((0, 0, 1))
    assert body.rotation_center.shape == (3,)
    np.testing.assert_allclose(body.rotation_center, np.array([0, 0, 1]))


def test_rotate_rotation_center_defined_as_tuple():
    body = cpt.FloatingBody()
    body.rotation_center = (1, 0, 0)
    body.rotate_z(np.pi)
    assert body.rotation_center.shape == (3,)
    np.testing.assert_allclose(body.rotation_center, np.array([-1, 0, 0]), atol=1e-10)


@pytest.mark.parametrize("z_center", [0, 2, -2])
@pytest.mark.parametrize("as_collection_of_meshes", [True, False])
def test_clipping_of_dofs(z_center, as_collection_of_meshes):
    """Check that clipping a body with a dof is the same as clipping the body ant then adding the dof."""
    full_sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(center=(0, 0, z_center), axial_symmetry=as_collection_of_meshes), name="sphere")
    axis = cpt.Axis(point=(1, 0, 0), vector=(1, 0, 0))

    full_sphere.add_rotation_dof(axis, name="test_dof")
    clipped_sphere = full_sphere.keep_immersed_part(free_surface=0.0, water_depth=np.inf, inplace=False)

    other_clipped_sphere = cpt.FloatingBody(mesh=clipped_sphere.mesh, name="other_sphere")
    other_clipped_sphere.add_rotation_dof(axis, name="test_dof")

    if clipped_sphere.mesh.nb_faces > 0:
        assert np.allclose(clipped_sphere.dofs['test_dof'], other_clipped_sphere.dofs['test_dof'])
    else:
        assert len(clipped_sphere.dofs['test_dof']) == 0


def test_cropping_body_with_manual_dof():
    # https://github.com/capytaine/capytaine/issues/204
    sphere = cpt.FloatingBody(cpt.mesh_sphere())
    sphere.dofs["Surge"] = [(1, 0, 0) for face in sphere.mesh.faces]
    sphere.keep_immersed_part()


def test_immersed_part():
    full_sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(), name="ball")
    immersed_sphere = full_sphere.immersed_part()
    assert immersed_sphere is not full_sphere
    assert immersed_sphere.mesh == full_sphere.mesh.immersed_part()
    assert immersed_sphere.mesh.axis_aligned_bbox[5] <= 0.0
    assert immersed_sphere.name == "ball"
    full_sphere.translate_x(2.0)
    new_immersed_sphere = full_sphere.immersed_part()
    assert new_immersed_sphere is not immersed_sphere
    assert not np.allclose(new_immersed_sphere.mesh.axis_aligned_bbox, immersed_sphere.mesh.axis_aligned_bbox)


def test_mincing():
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=10, radius=0.5))
    body = body.minced((4, 1, 1))
    assert len(body.mesh) == 2
    assert np.all(body.mesh[0].faces_centers[:, 0] < 0)
    assert isinstance(body.mesh[0][0], cpt.Mesh)
    body = body.minced((1, 2, 2))
    assert isinstance(body.mesh[0][0][0][0], cpt.Mesh)


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
    my_axis = cpt.Axis((0, 1, 0), point=(0,0,0))
    sphere.add_rotation_dof(axis=my_axis)
    sphere.keep_immersed_part()

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


def test_cluster_bodies():
    centers = np.array([[0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [11.0, 0.0, 0.0],
                        [10.0, 1.0, 0.0]])
    meshes = [cpt.mesh_sphere(center=c, name=str(c)) for c in centers]
    bodies = [cpt.FloatingBody(mesh=m) for m in meshes]
    joined_bodies = cpt.FloatingBody.join_bodies(*bodies)
    clustered_bodies = cpt.FloatingBody.cluster_bodies(*bodies)
    assert isinstance(clustered_bodies, cpt.FloatingBody)
    assert isinstance(clustered_bodies.mesh, cpt.CollectionOfMeshes)
    assert clustered_bodies.mesh.merged() == joined_bodies.mesh.merged()
    assert meshes[0] in clustered_bodies.mesh  # The first body is at the top level independently from the other three
