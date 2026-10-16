import pytest
import numpy as np

from capytaine.bodies.dofs import TranslationDof, RotationDof, DofOnSubmesh
import capytaine_test_helpers as helpers

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


ALL_DOFS = [
    TranslationDof((0, 0, 1)),
    RotationDof((0, 0, 0), (1, 0, 0)),
    RotationDof((0, 0, 0), (0, 1, 0)),
    RotationDof((0, 0, 0), (0, 0, 1)),
]

@pytest.mark.parametrize("dof", ALL_DOFS)
@pytest.mark.parametrize("other_dof", ALL_DOFS)
def test_second_order_motion_symmetry(dof, other_dof):
    points = RNG.uniform(size=(10, 3))
    d2u = dof.evaluate_second_order_motion_at_points(other_dof, points)
    d2u_ = other_dof.evaluate_second_order_motion_at_points(dof, points)
    assert np.allclose(d2u, d2u_)


def wamit_second_order_matrix(a):
    # Matrix H of eq. (15.49) of the WAMIT manual v7.6 (roll-pitch-yaw sequence)
    return 0.5 * np.array([
        [-(a[1]**2 + a[2]**2), 0.0, 0.0],
        [2*a[0]*a[1], -(a[0]**2 + a[2]**2), 0.0],
        [2*a[0]*a[2], 2*a[1]*a[2], -(a[0]**2 + a[1]**2)],
    ])


def test_second_order_rotation_same_convention_as_wamit():
    center = np.array([0.5, -1.0, 2.0])
    points = RNG.uniform(size=(10, 3))
    dofs = [RotationDof(center, d) for d in np.eye(3)]  # roll, pitch, yaw
    amplitudes = RNG.uniform(size=3)
    second_order = 0.5 * sum(
        amplitudes[k] * amplitudes[j] * dofs[k].evaluate_second_order_motion_at_points(dofs[j], points)
        for k in range(3) for j in range(3)
    )
    expected = (points - center) @ wamit_second_order_matrix(amplitudes).T
    assert np.allclose(second_order, expected)


@pytest.mark.parametrize("dof", [TranslationDof((0, 0, 1)), RotationDof((0, 0, 0), (1, 0, 0))])
def test_second_order_motion_with_translation_is_zero(dof):
    points = RNG.uniform(size=(10, 3))
    translation = TranslationDof((1, 0, 0))
    assert np.allclose(dof.evaluate_second_order_motion_at_points(translation, points), 0.0)
    assert np.allclose(translation.evaluate_second_order_motion_at_points(dof, points), 0.0)


def test_second_order_rotation_different_centers():
    with pytest.raises(ValueError):
        RotationDof((0, 0, 0), (1, 0, 0)).evaluate_second_order_motion_at_points(
            RotationDof((1, 0, 0), (0, 1, 0)), RNG.uniform(size=(3, 3)))


@pytest.mark.parametrize("faces_a", [slice(0, 5), range(2, 7), [0, 2, 4, 6], "mask"])
@pytest.mark.parametrize("faces_b", [slice(3, 8), [4, 5, 6, 7], None])
def test_second_order_motion_on_submeshes(faces_a, faces_b):
    mesh = helpers.small_sphere_mesh()
    if isinstance(faces_a, str):
        faces_a = RNG.uniform(size=mesh.nb_faces) > 0.5
    roll = RotationDof((0, 0, 0), (1, 0, 0))
    pitch = RotationDof((0, 0, 0), (0, 1, 0))
    full = roll.evaluate_second_order_motion(pitch, mesh)
    assert np.abs(full).max() > 0.0

    expected_mask = np.zeros(mesh.nb_faces, dtype=bool)
    expected_mask[faces_a] = True
    sub_a = DofOnSubmesh(roll, faces_a)
    if faces_b is None:  # only one of the two dofs is on a submesh
        other, result = pitch, sub_a.evaluate_second_order_motion(pitch, mesh)
        reverse = pitch.evaluate_second_order_motion(sub_a, mesh)
    else:
        mask_b = np.zeros(mesh.nb_faces, dtype=bool)
        mask_b[faces_b] = True
        expected_mask &= mask_b
        sub_b = DofOnSubmesh(pitch, faces_b)
        result = sub_a.evaluate_second_order_motion(sub_b, mesh)
        reverse = sub_b.evaluate_second_order_motion(sub_a, mesh)

    expected = np.where(expected_mask[:, None], full, 0.0)
    assert result.shape == (mesh.nb_faces, 3)
    assert np.allclose(result, expected)
    assert np.allclose(reverse, expected)


def _rotation_dofs_and_translation():
    center = (0.5, -0.3, 1.0)
    return [RotationDof(center, d) for d in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]] + [TranslationDof((1, 0, 0))]


def _finite_difference_gradient(f, points, eps=1e-3):
    # Jacobian convention: grad[..., i_dir, i_deriv_dir]
    grad = np.zeros((*points.shape[:-1], 3, 3))
    for j in range(3):
        step = np.zeros(3)
        step[j] = eps
        grad[..., :, j] = (f(points + step) - f(points - step)) / (2 * eps)
    return grad


@pytest.mark.parametrize("dof", _rotation_dofs_and_translation())
@pytest.mark.parametrize("other_dof", _rotation_dofs_and_translation())
def test_gradient_of_second_order_motion_finite_differences(dof, other_dof):
    points = RNG.uniform(size=(6, 3))
    grad = dof.evaluate_gradient_of_second_order_motion_at_points(other_dof, points)
    assert grad.shape == (6, 3, 3)
    expected = _finite_difference_gradient(lambda x: dof.evaluate_second_order_motion_at_points(other_dof, x), points)
    assert np.allclose(grad, expected, atol=1e-8)


def test_gradient_of_second_order_motion_is_not_trivial_and_symmetric():
    roll, pitch, yaw, _ = _rotation_dofs_and_translation()
    points = RNG.uniform(size=(2, 3))
    assert np.abs(roll.evaluate_gradient_of_second_order_motion_at_points(pitch, points)).max() > 0.0
    assert np.allclose(roll.evaluate_gradient_of_second_order_motion_at_points(pitch, points),
                       pitch.evaluate_gradient_of_second_order_motion_at_points(roll, points))


def test_gradient_of_second_order_rotation_different_centers():
    with pytest.raises(ValueError):
        RotationDof((0, 0, 0), (1, 0, 0)).evaluate_gradient_of_second_order_motion_at_points(
            RotationDof((1, 0, 0), (0, 1, 0)), RNG.uniform(size=(3, 3)))


@pytest.mark.parametrize("faces_a", [slice(0, 5), [0, 2, 4, 6]])
@pytest.mark.parametrize("faces_b", [slice(3, 8), None])
def test_gradient_of_second_order_motion_on_submeshes(faces_a, faces_b):
    mesh = helpers.small_sphere_mesh()
    roll = RotationDof((0, 0, 0), (1, 0, 0))
    pitch = RotationDof((0, 0, 0), (0, 1, 0))
    full = roll.evaluate_gradient_of_second_order_motion(pitch, mesh)
    assert full.shape == (mesh.nb_faces, 3, 3)
    assert np.abs(full).max() > 0.0
    expected_fd = _finite_difference_gradient(
        lambda x: roll.evaluate_second_order_motion_at_points(pitch, x), mesh.faces_centers)
    assert np.allclose(full, expected_fd, atol=1e-8)

    expected_mask = np.zeros(mesh.nb_faces, dtype=bool)
    expected_mask[faces_a] = True
    sub_a = DofOnSubmesh(roll, faces_a)
    if faces_b is None:
        result = sub_a.evaluate_gradient_of_second_order_motion(pitch, mesh)
        reverse = pitch.evaluate_gradient_of_second_order_motion(sub_a, mesh)
    else:
        mask_b = np.zeros(mesh.nb_faces, dtype=bool)
        mask_b[faces_b] = True
        expected_mask &= mask_b
        sub_b = DofOnSubmesh(pitch, faces_b)
        result = sub_a.evaluate_gradient_of_second_order_motion(sub_b, mesh)
        reverse = sub_b.evaluate_gradient_of_second_order_motion(sub_a, mesh)
    expected = np.where(expected_mask[:, None, None], full, 0.0)
    assert np.allclose(result, expected)
    assert np.allclose(reverse, expected)
