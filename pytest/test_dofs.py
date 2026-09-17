import pytest
import numpy as np

from capytaine.bodies.dofs import TranslationDof, RotationDof, DofOnSubmesh

RNG = np.random.default_rng()

@pytest.mark.parametrize("shape", [(5, 3), (5, 1, 3), (5, 1, 1, 3)])
@pytest.mark.parametrize("dof", [
    TranslationDof((0, 0, 1)),
    RotationDof((0, 0, 0), (1, 0, 0)),
])
def test_evaluate_higher_dim_arrays(dof, shape):
    points = RNG.uniform(size=shape)
    motions = dof.evaluate_motion_at_points(points)
    assert motions.shape == shape
    grad = dof.evaluate_gradient_of_motion_at_points(points)
    assert grad.shape == (*shape, 3)
