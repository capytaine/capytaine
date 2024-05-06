import pytest

import numpy as np
import capytaine as cpt


@pytest.fixture
def body_without_lid():
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    body_without_lid = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())
    return body_without_lid


@pytest.fixture
def body_with_lid():
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    lid_mesh = cpt.mesh_rectangle(size=(1, 1), resolution=(4, 4), center=(0, 0, 0), normal=(0, 0, 1))
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

    params = [0.0, np.inf, 1.0, solver.green_function]
    S_with, K_with = solver.engine.build_matrices(
            body_with_lid.mesh_including_lid, body_with_lid.mesh_including_lid,
            *params, adjoint_double_layer=True,
            )
    S_without, K_without = solver.engine.build_matrices(
            body_without_lid.mesh, body_without_lid.mesh,
            *params, adjoint_double_layer=True,
            )
    _, D_with = solver.engine.build_matrices(
            body_with_lid.mesh_including_lid, body_with_lid.mesh_including_lid,
            *params, adjoint_double_layer=False,
            )
    _, D_without = solver.engine.build_matrices(
            body_without_lid.mesh, body_without_lid.mesh,
            *params, adjoint_double_layer=False,
            )

    np.testing.assert_allclose(K_with[:n_hull_mesh, :n_hull_mesh], K_without, atol=1e-8)
    np.testing.assert_allclose(D_with[:n_hull_mesh, :n_hull_mesh], D_without, atol=1e-8)
    np.testing.assert_allclose(S_with[:n_hull_mesh, :n_hull_mesh], S_without, atol=1e-8)


def test_effect_of_lid_on_regular_frequency_diffraction_force(body_without_lid, body_with_lid):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))

    pb_with = cpt.DiffractionProblem(body=body_with_lid, wavelength=3.0)
    res_with = solver.solve(pb_with)
    f_with = res_with.forces["Heave"]

    pb_without = cpt.DiffractionProblem(body=body_without_lid, wavelength=3.0)
    res_without = solver.solve(pb_without)
    f_without = res_without.forces["Heave"]

    assert f_with == pytest.approx(f_without, rel=5e-2)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
def test_panel_on_free_surface(water_depth):
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.5))
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    body.mesh.compute_quadrature("Gauss-Legendre 2")
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, water_depth=water_depth, radiating_dof="Heave")
    res = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).solve(pb)
