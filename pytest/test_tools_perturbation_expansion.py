import numpy as np
import pytest

from capytaine.tools.perturbation_expansion import PerturbationExpansion, DifferenceFrequencyExpansion


def test_addition_scalar():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    result = a + b
    assert result.order == 2
    assert result.coefficients == (5, 7, 9)


def test_subtraction_and_negation():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    assert (a - b).coefficients == (-3, -3, -3)
    assert (-a).coefficients == (-1, -2, -3)
    assert (10 - a).coefficients == (9, -2, -3)


def test_elementwise_multiplication_truncated_cauchy_product():
    a = PerturbationExpansion([1, 2, 3])
    b = PerturbationExpansion([4, 5, 6])
    result = a * b
    assert result.order == 2
    assert result.coefficients == (4, 13, 28)


def test_scalar_broadcasting():
    a = PerturbationExpansion([np.array([1.0, 2.0]), np.array([3.0, 4.0])])
    result = a + 1.0
    np.testing.assert_allclose(result.coefficients[0], [2.0, 3.0])
    np.testing.assert_allclose(result.coefficients[1], [3.0, 4.0])

    result2 = 2.0 * a
    np.testing.assert_allclose(result2.coefficients[0], [2.0, 4.0])
    np.testing.assert_allclose(result2.coefficients[1], [6.0, 8.0])


def test_order_mismatch_raises():
    a = PerturbationExpansion([1, 2])
    b = PerturbationExpansion([1, 2, 3])
    with pytest.raises(ValueError):
        a + b
    with pytest.raises(ValueError):
        a * b


def test_getitem_indexes_into_contained_arrays():
    a = PerturbationExpansion([np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])])

    element = a[1]
    assert isinstance(element, PerturbationExpansion)
    assert element.order == a.order
    assert element.coefficients[0] == 2.0
    assert element.coefficients[1] == 5.0

    sliced = a[1:]
    np.testing.assert_allclose(sliced.coefficients[0], [2.0, 3.0])
    np.testing.assert_allclose(sliced.coefficients[1], [5.0, 6.0])

    with pytest.raises(IndexError):
        a[10]


def test_len_and_iter_act_on_contained_arrays():
    a = PerturbationExpansion([np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])])
    assert len(a) == 3

    elements = list(a)
    assert len(elements) == 3
    for i, element in enumerate(elements):
        assert isinstance(element, PerturbationExpansion)
        assert element.coefficients[0] == a.coefficients[0][i]
        assert element.coefficients[1] == a.coefficients[1][i]


def test_matmul_cauchy_product():
    I = np.eye(2)
    A1 = np.array([[0.0, 1.0], [2.0, 0.0]])
    A2 = np.array([[1.0, 0.0], [0.0, 1.0]])
    x = PerturbationExpansion([I, A1, A2])

    B0 = np.array([[2.0, 0.0], [0.0, 2.0]])
    B1 = np.array([[0.0, 0.0], [1.0, 0.0]])
    B2 = np.array([[0.0, 1.0], [0.0, 0.0]])
    y = PerturbationExpansion([B0, B1, B2])

    result = x @ y
    expected0 = I @ B0
    expected1 = I @ B1 + A1 @ B0
    expected2 = I @ B2 + A1 @ B1 + A2 @ B0
    np.testing.assert_allclose(result.coefficients[0], expected0)
    np.testing.assert_allclose(result.coefficients[1], expected1)
    np.testing.assert_allclose(result.coefficients[2], expected2)


def test_evaluate():
    a = PerturbationExpansion([1.0, 2.0, 3.0])
    assert a.evaluate(0.1) == pytest.approx(1.0 + 2.0 * 0.1 + 3.0 * 0.1**2)


def test_inverse_deformation_cofactor():
    I = np.eye(3)
    nabla_u = np.random.rand(3, 3)
    J = PerturbationExpansion([I, nabla_u])
    A = J.det() * J.inv().T
    expected1 = np.trace(nabla_u) * np.eye(3) - nabla_u.T
    assert np.allclose(A.coefficients[1], expected1)


def random_matrix_expansion(shape, order=2, seed=0):
    rng = np.random.default_rng(seed)
    coefficients = [np.broadcast_to(np.eye(3), shape + (3, 3))]
    coefficients += [rng.random(shape + (3, 3)) for _ in range(order)]
    return PerturbationExpansion(coefficients)


def test_radd_rmul_and_python_sum():
    a = PerturbationExpansion([np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([5.0, 6.0])])
    assert sum([a, a, a]).coefficients[2].tolist() == [15.0, 18.0]
    np.testing.assert_allclose((np.array([1.0, 10.0]) * a).coefficients[1], [3.0, 40.0])
    np.testing.assert_allclose((np.array([1.0, 10.0]) + a).coefficients[0], [2.0, 12.0])
    np.testing.assert_allclose((np.float64(2.0) * a).coefficients[2], [10.0, 12.0])
    np.testing.assert_allclose((np.array([1.0, 10.0]) - a).coefficients[1], [-3.0, -4.0])


def test_constant_matmul_on_both_sides():
    A = np.random.default_rng(1).random((3, 3))
    x = random_matrix_expansion(())
    left = A @ x
    right = x @ A
    for k in range(3):
        np.testing.assert_allclose(left.coefficients[k], A @ x.coefficients[k])
        np.testing.assert_allclose(right.coefficients[k], x.coefficients[k] @ A)


@pytest.mark.parametrize("shape", [(), (5,), (2, 5)])
def test_batched_det_and_inv(shape):
    F = random_matrix_expansion(shape)
    det, inv = F.det(), F.inv()
    assert det.order == 2 and inv.order == 2
    errors = []
    for eps in [1e-2, 5e-3]:
        M = F.evaluate(eps)
        errors.append((
            np.max(np.abs(det.evaluate(eps) - np.linalg.det(M))),
            np.max(np.abs(inv.evaluate(eps) - np.linalg.inv(M))),
        ))
    assert np.shape(det.coefficients[1]) == shape
    assert inv.coefficients[1].shape == shape + (3, 3)
    for i in range(2):
        assert errors[0][i] < 20 * 1e-2**3
        # Error is O(eps^3): halving eps divides it by about 8
        assert errors[0][i] / errors[1][i] == pytest.approx(8.0, rel=0.3)


def test_inv_and_det_at_order_2():
    F = random_matrix_expansion(())
    A1, A2 = F.coefficients[1], F.coefficients[2]
    np.testing.assert_allclose(F.inv().coefficients[2], -A2 + A1 @ A1)
    expected_det2 = np.trace(A2) + 0.5 * (np.trace(A1)**2 - np.trace(A1 @ A1))
    np.testing.assert_allclose(F.det().coefficients[2], expected_det2)


def test_det_inv_requires_identity_at_order_0():
    F = PerturbationExpansion([2 * np.eye(3), np.ones((3, 3))])
    with pytest.raises(NotImplementedError):
        F.det()
    with pytest.raises(NotImplementedError):
        F.inv()


def test_nanson_cofactor():
    F = random_matrix_expansion((4,))
    A = F.det().expand_dims(-1).expand_dims(-1) * F.inv().T
    nabla_u = F.coefficients[1]
    expected1 = np.trace(nabla_u, axis1=-2, axis2=-1)[:, None, None] * np.eye(3) - np.swapaxes(nabla_u, -1, -2)
    np.testing.assert_allclose(A.coefficients[1], expected1)


def test_field_alignment_needs_explicit_expand_dims():
    # A (nb_faces,) field and a (nb_faces, 3, 3) field: NumPy broadcasting pads in front and fails
    F = random_matrix_expansion((4,))
    d = F.det()
    with pytest.raises(ValueError):
        d * F.inv().T
    explicit = d.expand_dims(-1).expand_dims(-1) * F.inv().T
    assert explicit.coefficients[2].shape == (4, 3, 3)


def test_trace_sum_expand_dims_squeeze():
    rng = np.random.default_rng(2)
    c = [rng.random((4, 3, 3)) for _ in range(3)]
    x = PerturbationExpansion(c)
    for k in range(3):
        np.testing.assert_allclose(x.trace().coefficients[k], np.trace(c[k], axis1=-2, axis2=-1))
        np.testing.assert_allclose(x.sum(-3).coefficients[k], np.sum(c[k], axis=0))
        np.testing.assert_allclose(x.sum((-2, -1)).coefficients[k], np.sum(c[k], axis=(1, 2)))
        np.testing.assert_allclose(x.T.coefficients[k], np.swapaxes(c[k], -1, -2))
    e = x.expand_dims(-3)
    assert [co.shape for co in e.coefficients] == [(4, 1, 3, 3)] * 3
    assert e.squeeze(-3).coefficients[1].shape == (4, 3, 3)
    for method in [x.sum, x.expand_dims, x.squeeze]:
        with pytest.raises(ValueError, match="negative axis"):
            method(0)


def test_indexing_variants():
    rng = np.random.default_rng(3)
    c = [rng.random((4, 3)) for _ in range(3)]
    x = PerturbationExpansion(c)
    for index in [1, slice(1, 3), (slice(None), 2), (Ellipsis, 0), (1, 2), (slice(None), None)]:
        y = x[index]
        for k in range(3):
            np.testing.assert_allclose(y.coefficients[k], c[k][index])


def test_broadcast_dofs_constant_against_expansion():
    nb_dofs, nb_faces = 6, 4
    rng = np.random.default_rng(4)
    constant = rng.random((nb_dofs, nb_faces, 3))
    x = PerturbationExpansion([rng.random((nb_faces, 3)) for _ in range(3)])
    for result in [constant * x, x * constant, constant + x, x - constant, constant - x]:
        assert result.coefficients[0].shape == (nb_dofs, nb_faces, 3)
    product = constant * x
    for k in range(3):
        np.testing.assert_allclose(product.coefficients[k], constant * x.coefficients[k])
    assert (constant * x).sum(-1).coefficients[2].shape == (nb_dofs, nb_faces)


def test_matrix_vector_product_batched():
    rng = np.random.default_rng(5)
    M = PerturbationExpansion([rng.random((4, 3, 3)) for _ in range(3)])
    v = PerturbationExpansion([rng.random((4, 3)) for _ in range(3)])
    result = (M @ v.expand_dims(-1)).squeeze(-1)
    expected2 = sum(
        np.einsum("fab,fb->fa", M.coefficients[i], v.coefficients[2 - i]) for i in range(3)
    )
    np.testing.assert_allclose(result.coefficients[2], expected2)


# DifferenceFrequencyExpansion
# The reference is computed in the time domain: the real signals of the inputs are sampled,
# the formula is applied to them with a PerturbationExpansion, and the Fourier coefficients
# of its order 2 are extracted at the difference frequencies omega_m - omega_n.

OMEGAS = np.array([10.0, 11.0, 13.0])  # No difference frequency equals another one or a sum frequency
TIMES = np.linspace(0.0, 2*np.pi, 128, endpoint=False)  # A common period of all the frequencies


def random_wave_expansion(field_shape, seed, order_0=None, order_2=False):
    rng = np.random.default_rng(seed)
    nb_waves = len(OMEGAS)
    x0 = rng.random(field_shape) if order_0 is None else order_0
    x1 = rng.random((nb_waves,) + field_shape) + 1j*rng.random((nb_waves,) + field_shape)
    x2 = 0.0
    if order_2:
        q = rng.random((nb_waves, nb_waves) + field_shape) + 1j*rng.random((nb_waves, nb_waves) + field_shape)
        x2 = (q + np.conj(np.swapaxes(q, 0, 1))) / 2  # Hermitian
    return DifferenceFrequencyExpansion([x0, x1, x2])


def time_signal(x, t):
    """Real PerturbationExpansion of the signal represented by `x`, at time `t`."""
    phases = np.exp(-1j*OMEGAS*t)
    x1 = np.real(np.tensordot(phases, x.coefficients[1], axes=1))
    x2 = np.real(np.tensordot(np.outer(phases, np.conj(phases)), x.coefficients[2], axes=2))
    return PerturbationExpansion([x.coefficients[0], x1, x2])


def time_domain_qtf(f, *inputs):
    samples = np.array([
        f(*(time_signal(x, t) if isinstance(x, DifferenceFrequencyExpansion) else x for x in inputs)).coefficients[2]
        for t in TIMES
    ])
    nb_waves = len(OMEGAS)
    qtf = np.empty((nb_waves, nb_waves) + samples.shape[1:], dtype=complex)
    for m in range(nb_waves):
        for n in range(nb_waves):
            phases = np.exp(1j*(OMEGAS[m] - OMEGAS[n])*TIMES).reshape((-1,) + (1,)*(samples.ndim - 1))
            qtf[m, n] = np.mean(samples * phases, axis=0)
    return qtf


def check_against_time_domain(f, *inputs):
    result = f(*inputs)
    assert isinstance(result, DifferenceFrequencyExpansion)
    qtf = result.coefficients[2]
    reference = time_domain_qtf(f, *inputs)
    assert qtf.shape == reference.shape
    off_diagonal = ~np.eye(len(OMEGAS), dtype=bool)
    np.testing.assert_allclose(qtf[off_diagonal], reference[off_diagonal], rtol=1e-10, atol=1e-10)
    # The diagonal terms all are at frequency 0: only their sum (the mean) can be extracted.
    np.testing.assert_allclose(np.trace(qtf), reference[0, 0], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(qtf, np.conj(np.swapaxes(qtf, 0, 1)), atol=1e-12)


def test_difference_frequency_product():
    a = random_wave_expansion((4,), seed=0)
    b = random_wave_expansion((4,), seed=1, order_2=True)
    check_against_time_domain(lambda a, b: a * b, a, b)
    check_against_time_domain(lambda a, b: (a - 2.0) * (b + 1.0) * a, a, b)


def test_difference_frequency_diagonal_is_the_mean():
    a = random_wave_expansion((4,), seed=0, order_0=0.0)
    b = random_wave_expansion((4,), seed=1, order_0=0.0)
    qtf = (a * b).coefficients[2]
    for m in range(len(OMEGAS)):
        np.testing.assert_allclose(qtf[m, m], 0.5 * np.real(a.coefficients[1][m] * np.conj(b.coefficients[1][m])))


def test_difference_frequency_matrix_product():
    A = random_wave_expansion((4, 3, 3), seed=0)
    B = random_wave_expansion((4, 3, 3), seed=1, order_2=True)
    check_against_time_domain(lambda A, B: A @ B, A, B)
    check_against_time_domain(lambda A, B: A.T @ B @ np.ones((3, 3)), A, B)


def test_difference_frequency_det_inv():
    F = random_wave_expansion((4, 3, 3), seed=0, order_0=np.eye(3), order_2=True)
    check_against_time_domain(lambda F: F.det(), F)
    check_against_time_domain(lambda F: F.inv(), F)
    check_against_time_domain(lambda F: F.det().expand_dims(-1).expand_dims(-1) * F.inv().T, F)


def test_difference_frequency_with_constants_of_other_shapes():
    nb_dofs, nb_faces = 2, 4
    rng = np.random.default_rng(6)
    motion = rng.random((nb_dofs, nb_faces, 3))  # Constant with more field axes than the expansions
    gradient_of_motion = rng.random((nb_dofs, nb_faces, 3, 3))
    u = random_wave_expansion((nb_faces, 3), seed=0, order_2=True)
    p = random_wave_expansion((nb_faces,), seed=1)
    q = random_wave_expansion((), seed=2, order_0=0.0)  # Scalar field
    check_against_time_domain(
        lambda u, p: (p * ((motion + (gradient_of_motion @ u.expand_dims(-1)).squeeze(-1)) * u).sum(-1)).sum(-1),
        u, p)
    check_against_time_domain(lambda q, u: q * motion + u, q, u)
    assert (q * motion + u).coefficients[2].shape == (3, 3, nb_dofs, nb_faces, 3)


def test_difference_frequency_shapes_and_indexing():
    x = DifferenceFrequencyExpansion([0.0, np.ones((3, 4, 2)), 0.0])
    assert x.nb_waves == 3
    assert [c.shape for c in x.coefficients] == [(4, 2), (3, 4, 2), (3, 3, 4, 2)]
    assert len(x) == 4
    assert [c.shape for c in x[..., 1].coefficients] == [(4,), (3, 4), (3, 3, 4)]
    assert [c.shape for c in x.sum(-2).coefficients] == [(2,), (3, 2), (3, 3, 2)]
    with pytest.raises(IndexError):
        x[1]


def test_difference_frequency_errors():
    x = random_wave_expansion((4, 3, 3), seed=0)
    with pytest.raises(ValueError, match="real"):
        DifferenceFrequencyExpansion([1j*np.ones(3), np.ones((2, 3)), 0.0])
    with pytest.raises(ValueError, match="real"):
        1j * x  # Complex factors apply to the first-order amplitudes only
    with pytest.raises(ValueError, match="waves"):
        x + DifferenceFrequencyExpansion([0.0, np.ones((2, 4, 3, 3)), 0.0])
    with pytest.raises(TypeError):
        x + PerturbationExpansion([np.eye(3), np.eye(3), np.eye(3)])
    with pytest.raises(TypeError):
        PerturbationExpansion([np.eye(3), np.eye(3), np.eye(3)]) + x
    v = random_wave_expansion((3,), seed=1)  # Vector field
    with pytest.raises(ValueError, match="2 field axes"):
        x @ v
    with pytest.raises(ValueError, match="2 field axes"):
        np.ones(3) @ x
    with pytest.raises(NotImplementedError):
        x.evaluate(0.1)
