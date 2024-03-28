import pytest

import numpy as np
import capytaine as cpt
from capytaine.bem.airy_waves import froude_krylov_force


def test_Froude_Krylov():
    sphere = cpt.FloatingBody(cpt.mesh_sphere(radius=1.0, resolution=(6, 12)).immersed_part())
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, water_depth=np.inf)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 27596, rtol=1e-3)

    problem = cpt.DiffractionProblem(body=sphere, omega=2.0, water_depth=np.inf)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 22491, rtol=1e-3)

    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, water_depth=10.0)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 27610, rtol=1e-3)
