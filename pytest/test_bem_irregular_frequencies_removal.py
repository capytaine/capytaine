# import pytest

import numpy as np
import capytaine as cpt


def test_panel_on_free_surface():
    mesh = cpt.mesh_parallelepiped(center=(0.0, 0.0, -0.5))
    body = cpt.FloatingBody(mesh, cpt.rigid_body_dofs())
    body.mesh.compute_quadrature("Gauss-Legendre 2")
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    res = cpt.BEMSolver().solve(pb)
