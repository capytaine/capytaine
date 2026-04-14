import pytest

import numpy as np
import capytaine as cpt
import xarray as xr

from capytaine.bodies.multibodies import Multibody


def test_multibody_without_names():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            name="body"
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            name="body"
            )
    with pytest.raises(ValueError):
        Multibody([body_1, body_2],)


def test_multibody_resolution():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            name="body_2",
            )
    multi = Multibody(
        [body_1, body_2],
        # own_dofs={
        #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
        # }
    )
    solver = cpt.BEMSolver()
    res = solver.solve(cpt.DiffractionProblem(body=multi, omega=1.0, wave_direction=0.0))
    ref_res = solver.solve(cpt.DiffractionProblem(body=multi.as_FloatingBody(), omega=1.0, wave_direction=0.0))
    assert all(np.isclose(ref_res.forces[k], res.forces[k]) for k in multi.dofs)

    res = solver.solve(cpt.RadiationProblem(body=multi, omega=1.0, radiating_dof="body_1__Heave"))
    ref_res = solver.solve(cpt.RadiationProblem(body=multi.as_FloatingBody(), omega=1.0, radiating_dof="body_1__Heave"))
    assert all(np.isclose(ref_res.forces[k], res.forces[k]) for k in multi.dofs)

def test_multibody_hydrostatics():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0),
            mass=1000.0,
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            center_of_mass=(2, 0, 0),
            mass=1000.0,
            name="body_2",
            )
    multi = Multibody(
            [body_1, body_2],
            # own_dofs={
            #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
            #     }
            )
    assert np.allclose(multi.center_of_mass['body_1'], body_1.center_of_mass)
    assert np.allclose(multi.center_of_mass['body_2'], body_2.center_of_mass)
    assert np.allclose(multi.center_of_buoyancy['body_1'], body_1.center_of_buoyancy)
    assert np.allclose(multi.center_of_buoyancy['body_2'], body_2.center_of_buoyancy)

def test_hydrostatic_stiffness():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_parallelepiped(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0),
            mass=1000.0,
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_parallelepiped(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            center_of_mass=(2, 0, 0),
            mass=1000.0,
            name="body_2",
            )
    multi = Multibody(
            [body_1, body_2],
            # own_dofs={
            #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
            #     }
            )
    h = multi.compute_hydrostatic_stiffness()
    assert np.allclose(h.values[:6, :6], body_1.compute_hydrostatic_stiffness().values)
    assert np.allclose(h.values[6:, 6:], body_2.compute_hydrostatic_stiffness().values)
    assert np.allclose(h.values[6:, :6], 0.0)
    assert np.allclose(h.values[:6, 6:], 0.0)

    m = multi.compute_rigid_body_inertia()
    assert np.allclose(m.values[:6, :6], body_1.compute_rigid_body_inertia().values)
    assert np.allclose(m.values[6:, 6:], body_2.compute_rigid_body_inertia().values)
    assert np.allclose(m.values[6:, :6], 0.0)
    assert np.allclose(m.values[:6, 6:], 0.0)


def test_assemble_regular_array():
    body = cpt.FloatingBody(mesh=cpt.mesh_sphere())
    body.add_all_rigid_body_dofs()
    array = body.assemble_regular_array(distance=2.0, nb_bodies=(2, 3))
    assert array.mesh.nb_faces == 6*body.mesh.nb_faces

    # Check name and order of the dofs
    assert list(array.dofs.keys())[0:3] == ["0_0__Surge", "0_0__Sway", "0_0__Heave"]
    assert "2_1__Heave" not in array.dofs.keys()

    # Check that the dofs corresponds to the right panels
    faces_1_0 = np.where(array.dofs["1_0__Heave"].evaluate_motion(array.mesh) != 0.0)[0]
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
    for dof_name, dof in fb_array.dofs.items():
        num_active_faces.append(np.count_nonzero(np.count_nonzero(dof.evaluate_motion(fb_array.mesh), axis=1)))

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


def test_clip_multibody():
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
    heave1 = both.dofs["body_1__Heave"].evaluate_motion(both.mesh)
    heave2 = both.dofs["body_2__Heave"].evaluate_motion(both.mesh)
    # heave1 and heave2 should be non-zero on half the faces
    assert np.count_nonzero(heave1[:, 2]) == both.mesh.nb_faces // 2
    assert np.count_nonzero(heave2[:, 2]) == both.mesh.nb_faces // 2

    immersed_both = both.immersed_part()
    heave1 = immersed_both.dofs["body_1__Heave"].evaluate_motion(immersed_both.mesh)
    heave2 = immersed_both.dofs["body_2__Heave"].evaluate_motion(immersed_both.mesh)
    # heave1 and heave2 should be non-zero on half the faces
    assert np.count_nonzero(heave1[:, 2]) == immersed_both.mesh.nb_faces // 2
    assert np.count_nonzero(heave2[:, 2]) == immersed_both.mesh.nb_faces // 2


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
    assert both.dofs["body_1__Heave"].evaluate_motion(both.mesh).shape[0] == both.mesh.nb_faces


@pytest.mark.parametrize("sym_2", [True, False])
@pytest.mark.parametrize("sym_1", [True, False])
@pytest.mark.parametrize("lid_2", [True, False])
@pytest.mark.parametrize("lid_1", [True, False])
def test_with_and_without_symmetry_with_and_without_lid(lid_1, lid_2, sym_1, sym_2):
    mesh_1 = cpt.mesh_horizontal_cylinder(reflection_symmetry=sym_1).immersed_part()
    body_1 = cpt.FloatingBody(
        mesh_1,
        lid_mesh=mesh_1.generate_lid(faces_max_radius=0.3) if lid_1 else None,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        name="body_1"
    )
    assert body_1.hull_mask.shape == (body_1.mesh_including_lid.nb_faces,)
    mesh_2 = cpt.mesh_horizontal_cylinder(center=(12, 0, 0), reflection_symmetry=sym_2).immersed_part()
    body_2 = cpt.FloatingBody(
        mesh_2,
        lid_mesh=mesh_2.generate_lid(faces_max_radius=0.3) if lid_2 else None,
        dofs=cpt.rigid_body_dofs(rotation_center=(12, 0, 0)),
        name="body_2"
    )
    assert body_2.hull_mask.shape == (body_2.mesh_including_lid.nb_faces,)
    both = cpt.Multibody([body_1, body_2])
    assert both.mesh.nb_faces == body_1.mesh.nb_faces + body_2.mesh.nb_faces
    def nb_faces(lid_mesh):
        if lid_mesh is None:
            return 0
        else:
            return lid_mesh.nb_faces
    assert nb_faces(both.lid_mesh) == nb_faces(body_1.lid_mesh) + nb_faces(body_2.lid_mesh)
    assert both.mesh_including_lid.nb_faces == body_1.mesh_including_lid.nb_faces + body_2.mesh_including_lid.nb_faces
    assert both.hull_mask.shape[0] == body_1.hull_mask.shape[0] + body_2.hull_mask.shape[0]
    assert np.allclose(both.mesh_including_lid.faces_centers[~both.hull_mask, 2], 0.0)
    pb = cpt.RadiationProblem(body=both, omega=1.0, radiating_dof='body_1__Surge')
    solver = cpt.BEMSolver()
    res = solver.solve(pb)
    both.dofs['body_1__Surge']
    assert res.forces['body_1__Surge'] == pytest.approx(2151.7+54.053j, rel=0.05)
