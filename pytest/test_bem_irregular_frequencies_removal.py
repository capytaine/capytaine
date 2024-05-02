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


def test_lid_mesh(body_with_lid):
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))
    pb = cpt.DiffractionProblem(body=body_with_lid, wavelength=2.0)
    solver.solve(pb)
    pb = cpt.RadiationProblem(body=body_with_lid, wavelength=2.0)
    solver.solve(pb)


def test_froude_krylov_force(body_without_lid, body_with_lid):
    # Froude-Krylov force should be unchanged by the lid
    from capytaine.bem.airy_waves import froude_krylov_force
    pb_without_lid = cpt.DiffractionProblem(body=body_without_lid)
    pb_with_lid = cpt.DiffractionProblem(body=body_with_lid)

    fk_with_lid = froude_krylov_force(pb_with_lid)["Heave"]
    fk_without_lid = froude_krylov_force(pb_without_lid)["Heave"]
    assert fk_with_lid == pytest.approx(fk_without_lid)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
def test_panel_on_free_surface(water_depth):
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.5))
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    body.mesh.compute_quadrature("Gauss-Legendre 2")
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, water_depth=water_depth, radiating_dof="Heave")
    res = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).solve(pb)
