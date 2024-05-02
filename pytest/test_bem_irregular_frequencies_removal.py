import pytest

import numpy as np
import capytaine as cpt


def test_lid_mesh():
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.4)).immersed_part()
    lid_mesh = cpt.mesh_rectangle(size=(1, 1), resolution=(4, 4), center=(0, 0, 0), normal=(0, 0, 1))
    body = cpt.FloatingBody(mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities='low_freq'))
    pb = cpt.DiffractionProblem(body=body, wavelength=2.0)
    solver.solve(pb)
    pb = cpt.RadiationProblem(body=body, wavelength=2.0)
    solver.solve(pb)


@pytest.mark.parametrize("water_depth", [np.inf, 10.0])
def test_panel_on_free_surface(water_depth):
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.5))
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    body.mesh.compute_quadrature("Gauss-Legendre 2")
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, water_depth=water_depth, radiating_dof="Heave")
    res = cpt.BEMSolver(green_function=cpt.Delhommeau(gf_singularities="low_freq")).solve(pb)
