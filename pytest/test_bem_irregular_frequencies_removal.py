import pytest

import numpy as np
import capytaine as cpt

from capytaine.new_meshes.symmetric_meshes import ReflectionSymmetricMesh
from capytaine.new_meshes.predefined import mesh_parallelepiped, mesh_rectangle, mesh_vertical_cylinder, mesh_horizontal_cylinder

def test_irr_freq_warning_no_lid():
    mesh = mesh_parallelepiped(size=(1, 1, 1)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=None, dofs=cpt.rigid_body_dofs())
    assert 6.0 < body.first_irregular_frequency_estimate() < 7.0


def test_irr_freq_warning_subsurface_lid():
    mesh = mesh_parallelepiped(size=(1, 1, 1)).immersed_part()
    lid_mesh = mesh_rectangle(size=(1, 1), center=(0, 0, -0.1))
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    assert 10.0 < body.first_irregular_frequency_estimate() < 11.0


def test_irr_freq_warning_surface_lid():
    mesh = mesh_parallelepiped(size=(1, 1, 1)).immersed_part()
    lid_mesh = mesh_rectangle(size=(1, 1), center=(0, 0, 0))
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    assert body.first_irregular_frequency_estimate() == np.inf


def test_lid_with_upwards_normals():
    mesh = mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    lid_mesh = mesh_rectangle(size=(1, 1), resolution=(4, 4), center=(0, 0, -0.1), normal=(0, 0, 1))
    # Like body_with_lid() but with upward normals
    body = cpt.FloatingBody(mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    np.testing.assert_allclose(body.lid_mesh.faces_normals[:, 2], -np.ones((body.lid_mesh.nb_faces,)))

    # Old lid_mesh is unchanged
    assert body.lid_mesh is not lid_mesh
    np.testing.assert_allclose(lid_mesh.faces_normals[:, 2], np.ones((body.lid_mesh.nb_faces,)))


@pytest.fixture
def body_without_lid():
    mesh = mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    body_without_lid = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())
    return body_without_lid


@pytest.fixture
def body_with_lid():
    mesh = mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    lid_mesh = mesh_rectangle(size=(1, 1), resolution=(4, 4), center=(0, 0, 0), normal=(0, 0, -1))
    body_with_lid = cpt.FloatingBody(mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    return body_with_lid


def test_effect_of_lid_on_froude_krylov_force(body_without_lid, body_with_lid):
    # Froude-Krylov force should be unchanged by the lid
    from capytaine.bem.airy_waves import froude_krylov_force
    pb_without_lid = cpt.DiffractionProblem(body=body_without_lid)
    pb_with_lid = cpt.DiffractionProblem(body=body_with_lid)

    fk_with_lid = froude_krylov_force(pb_with_lid)["Heave"]
    fk_without_lid = froude_krylov_force(pb_without_lid)["Heave"]
    assert fk_with_lid == pytest.approx(fk_without_lid)


def test_effect_of_lid_on_matrices(body_without_lid, body_with_lid):
    n_hull_mesh = body_without_lid.mesh.nb_faces
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))

    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
    S_with, K_with = solver.engine.build_matrices(
            body_with_lid.mesh_including_lid, body_with_lid.mesh_including_lid,
            **params, adjoint_double_layer=True,
            )
    S_without, K_without = solver.engine.build_matrices(
            body_without_lid.mesh, body_without_lid.mesh,
            **params, adjoint_double_layer=True,
            )
    _, D_with = solver.engine.build_matrices(
            body_with_lid.mesh_including_lid, body_with_lid.mesh_including_lid,
            **params, adjoint_double_layer=False,
            )
    _, D_without = solver.engine.build_matrices(
            body_without_lid.mesh, body_without_lid.mesh,
            **params, adjoint_double_layer=False,
            )

    np.testing.assert_allclose(K_with[:n_hull_mesh, :n_hull_mesh], K_without, atol=1e-8)
    np.testing.assert_allclose(D_with[:n_hull_mesh, :n_hull_mesh], D_without, atol=1e-8)
    np.testing.assert_allclose(S_with[:n_hull_mesh, :n_hull_mesh], S_without, atol=1e-8)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
@pytest.mark.parametrize("forward_speed", [0.0, 1.0])
def test_effect_of_lid_on_regular_frequency_diffraction_force(
        body_without_lid, body_with_lid, water_depth, forward_speed,
        ):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))

    pb_with = cpt.DiffractionProblem(
            body=body_with_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed,
            )
    res_with = solver.solve(pb_with)
    f_with = res_with.forces["Heave"]

    pb_without = cpt.DiffractionProblem(
            body=body_without_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed
            )
    res_without = solver.solve(pb_without)
    f_without = res_without.forces["Heave"]

    assert f_with == pytest.approx(f_without, rel=5e-2)


def test_effect_of_lid_on_infinite_frequency(
        body_without_lid, body_with_lid,
        ):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))

    pb_with = cpt.RadiationProblem(
            body=body_with_lid, wavenumber=np.inf, radiating_dof="Heave"
            )
    res_with = solver.solve(pb_with)
    f_with = res_with.forces["Heave"]

    pb_without = cpt.RadiationProblem(
            body=body_without_lid, wavenumber=np.inf, radiating_dof="Heave"
            )
    res_without = solver.solve(pb_without)
    f_without = res_without.forces["Heave"]

    assert f_with == pytest.approx(f_without, rel=5e-2)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
@pytest.mark.parametrize("forward_speed", [0.0, 1.0])
def test_effect_of_lid_on_regular_frequency_free_surface_elevation(
        body_without_lid, body_with_lid,
        water_depth, forward_speed,
        ):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))
    point = np.array([[4.0, 4.0]])

    pb_with = cpt.DiffractionProblem(
            body=body_with_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed,
            )
    res_with = solver.solve(pb_with)
    fse_with = solver.compute_free_surface_elevation(point, res_with)

    pb_without = cpt.DiffractionProblem(
            body=body_without_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed,
            )
    res_without = solver.solve(pb_without)
    fse_without = solver.compute_free_surface_elevation(point, res_without)

    assert fse_with == pytest.approx(fse_without, rel=5e-2)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
@pytest.mark.parametrize("forward_speed", [0.0, 1.0])
def test_effect_of_lid_on_regular_frequency_field_velocity(
        body_without_lid, body_with_lid,
        water_depth, forward_speed,
        ):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))
    point = np.array([[4.0, 4.0, -2.0]])

    pb_with = cpt.DiffractionProblem(
            body=body_with_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed,
            )
    res_with = solver.solve(pb_with)
    u_with = solver.compute_velocity(point, res_with)

    pb_without = cpt.DiffractionProblem(
            body=body_without_lid, wavelength=3.0,
            water_depth=water_depth, forward_speed=forward_speed,
            )
    res_without = solver.solve(pb_without)
    u_without = solver.compute_velocity(point, res_without)

    assert u_with == pytest.approx(u_without, rel=5e-2)


@pytest.mark.xfail  # Need to update body class
def test_lid_multibody(body_with_lid):
    two_bodies = body_with_lid + body_with_lid.translated_x(5.0)
    n = body_with_lid.mesh.nb_faces
    nl = body_with_lid.lid_mesh.nb_faces

    # Check that the lid are still positioned on top of each component
    def horizontal_center(mesh):
        return mesh.faces_centers.mean(axis=0)[:2]
    np.testing.assert_allclose(
        horizontal_center(two_bodies.lid_mesh.extract_faces(list(range(0, nl)))),
        horizontal_center(two_bodies.mesh.extract_faces(list(range(0, n)))),
        atol=1e-6
    )
    np.testing.assert_allclose(
        horizontal_center(two_bodies.lid_mesh.extract_faces(list(range(nl, 2*nl)))),
        horizontal_center(two_bodies.mesh.extract_faces(list(range(n, 2*n)))),
        atol=1e-6
    )

    pb = cpt.DiffractionProblem(body=two_bodies, wavelength=3.0)
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))
    solver.solve(pb)


def test_lid_with_plane_symmetry():
    mesh = mesh_horizontal_cylinder(reflection_symmetry=True).immersed_part()
    lid_mesh = ReflectionSymmetricMesh(
            mesh_rectangle(size=(1.0, 10,), faces_max_radius=0.5, center=(0, 0.5, -0.05,)),
            plane=mesh.plane
            )
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = cpt.BEMSolver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    from capytaine.tools.block_circulant_matrices import BlockCirculantMatrix
    assert isinstance(S, BlockCirculantMatrix)
    assert isinstance(K, BlockCirculantMatrix)

    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])


def test_lid_with_nested_plane_symmetry():
    mesh = mesh_parallelepiped(center=(0, 0, -1.0), reflection_symmetry=True)
    hull_mesh, lid_mesh = mesh.extract_lid()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = cpt.BEMSolver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    from capytaine.tools.block_circulant_matrices import BlockCirculantMatrix
    assert isinstance(S, BlockCirculantMatrix)
    assert isinstance(S.blocks[0], BlockCirculantMatrix)
    assert isinstance(K, BlockCirculantMatrix)
    assert isinstance(K.blocks[0], BlockCirculantMatrix)

    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])


def test_lid_with_rotation_symmetry():
    mesh = mesh_vertical_cylinder(length=2.0, center=(0, 0, -1.0), axial_symmetry=True)
    hull_mesh, lid_mesh = mesh.extract_lid()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = cpt.BEMSolver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    from capytaine.tools.block_circulant_matrices import BlockCirculantMatrix
    assert isinstance(S, BlockCirculantMatrix)
    assert isinstance(K, BlockCirculantMatrix)
    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
def test_panel_on_free_surface(water_depth):
    mesh = mesh_parallelepiped(center=(0.0, 0.0, -0.5)).with_quadrature("Gauss-Legendre 2")
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, water_depth=water_depth, radiating_dof="Heave")
    cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).solve(pb)


def test_panel_on_free_surface_with_high_freq(caplog):
    mesh = mesh_rectangle(center=(0.0, 0.0, 0.0), size=(1.0, 1.0))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).solve(pb)
    with caplog.at_level("WARNING"):
        cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="high_freq")).solve(pb)
    assert "free surface panel" in caplog.text
