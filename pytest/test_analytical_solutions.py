import pytest
import numpy as np

import capytaine as cpt
from capytaine.bem.airy_waves import froude_krylov_force

from scipy.special import j0, y0, jn, yn  # Bessel functions


@pytest.mark.parametrize("k", [5.0, 10.0])
def test_mccamy_and_fuchs(k):
    h = 1.0  # Water depth
    R = 0.3  # Cylinder radius

    mesh = cpt.mesh_vertical_cylinder(radius=R, length=h, center=(0, 0, -h/2), resolution=(0, 30, 20)).immersed_part(water_depth=h)
    body = cpt.FloatingBody(mesh=mesh, dofs={"Surge": np.array([[1.0, 0.0, 0.0] for face in mesh.faces])})
    pb = cpt.DiffractionProblem(body=body, wavenumber=k, wave_direction=0.0, water_depth=h)

    # Analytical solution of McCamy and Fuchs for wavenumber k
    m0 = pb.omega**2/pb.g
    dy = (y0(m0*R) - yn(2, m0*R))
    dj = (j0(m0*R) - jn(2, m0*R))
    analytical_force = (4*pb.rho*pb.g*h*R) * (-1j*2*h/R) * np.tanh(m0*h)/(m0*h)**2 * (dy + 1j*dj)/(dy**2 + dj**2)

    solver = cpt.BEMSolver()
    res = solver.solve(pb)
    numerical_force = res.force["Surge"] + froude_krylov_force(res)["Surge"]

    np.testing.assert_allclose(numerical_force, analytical_force, rtol=0.05)
