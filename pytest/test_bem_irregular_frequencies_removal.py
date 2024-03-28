import pytest

import numpy as np
import capytaine as cpt


@pytest.mark.parametrize("water_depth", [np.infty, 10.0])
def test_panel_on_free_surface(water_depth):
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.5))
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    body.mesh.compute_quadrature("Gauss-Legendre 2")
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, water_depth=water_depth, radiating_dof="Heave")
    res = cpt.BEMSolver().solve(pb)
